#!/bin/python
"""
    Male Relatedness Analyzer
    Calculates the relatedness through males, females and both in a specified group.
    Copyright (C) 2015  Nikhil Tej Kurapati
    nikhiltk@umich.edu

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import csv
import copy
import math
import itertools
import os
import sys
# remove 0relatedness from wegs and lineages
# change relatedness through mothers vs fathers


# Default row names
field_definitions = {"id": "UniqueName",
                     "father_id": "UniqueFather",
                     "mother_id": "UniqueMother",
                     "is_alive_str": "Alive",
                     "sex_str": "Sex3",
                     "birth_year_str": "YOB_combo",
                     "simple_name": "New_name",
                     "primary_group": "Pedigree.name"
                     }

try:
    with open("excluded_individuals.txt") as f:
        excludes = f.readlines()
except IOError:
    excludes = {}


def main():
    from argparse import ArgumentParser
    import os.path
    # Parsing arguments passed from the command line
    parser = ArgumentParser(description="Calculates intragroup and intergroup relatedness "
                                        "for males in groups from a spreadsheet of data.  "
                                        "Assumes data is a list of rows of males")
    parser.add_argument("-i", dest="filename", required=True, help="input spreadsheet in csv format")

    parser.add_argument("-a", dest="age_restriction", required=False, default=0,
                        help="Specify only analyzing individuals who were born after the year specified")
    parser.add_argument("-s", dest="sex_str", required=False, default="Sex3",
                        help="The column name specifying the sex of the individuals: 1=male 0=female")
    parser.add_argument("-uid", dest="id", required=False, default="UniqueName",
                        help="The column name specifying the individual's unique id")
    parser.add_argument("-fid", dest="father_id", required=False, default="UniqueFather",
                        help="The column name specifying the unique id of the father")
    parser.add_argument("-mid", dest="mother_id", required=False, default="UniqueMother",
                        help="The column name specifying the unique id of the mother")
    parser.add_argument("-y", dest="birth_year_str", required=False, default="YOB_combo",
                        help="The column name specifying the year of birth")
    parser.add_argument("-f", dest="filter_variable", required=False, default="",
                        help="The column name specifying which rows to use or not to use")
    parser.add_argument("-n", dest="name", required=False, default="New_name",
                        help="The column name specifying the name (not necessarily unique of the individual)")
    parser.add_argument("-p", dest="primary_group", required=True, default="",
                        help="The column name specifying the primary group of the individual")
    parser.add_argument("-sg", dest="secondary_group", required=False, default="",
                        help="The column name specifying the secondary group of the individual")
    parser.add_argument("-m", dest="fill_missing_mothers", default=0.0, type=float,
                        help="Assume paternal sibs have same mother?")
    parser.add_argument("-gf", dest="regex_group_filter", default="",
                        help="regex filter for primary group names")
    parser.add_argument("-g2f", dest="regex_group2_filter", default="",
                        help="regex filter for secondary group  names")
    parser.add_argument("-an", dest="analysis_name", default="",
                        help="name for this specific set of conditions to analyze")
    args = parser.parse_args()

    # Setting up arguments
    global analysis_name
    analysis_name = args.analysis_name
    global age_restriction
    age_restriction = args.age_restriction
    global filter_variable
    filter_variable = args.filter_variable
    global field_definitions

    global missing_mother_fill_rate
    global group_restriction
    group_restriction = args.regex_group_filter

    global group2_restriction
    group2_restriction = args.regex_group2_filter

    missing_mother_fill_rate = args.fill_missing_mothers
    field_definitions["sex_str"] = args.sex_str
    field_definitions["id"] = args.id
    field_definitions["father_id"] = args.father_id
    field_definitions["mother_id"] = args.mother_id
    field_definitions['birth_year_str'] = args.birth_year_str
    field_definitions["simple_name"] = args.name
    field_definitions['primary_group'] = args.primary_group
    field_definitions['secondary_group'] = args.secondary_group

    # Setting up output path
    output_dir = args.filename + '_data'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Do the actual analysis
    run_analysis(args.filename, output_dir, grouping_criteria=field_definitions['primary_group'])


def run_analysis(filename, output_dir, grouping_criteria='group'):
    """
    Runs the actual analysis
    :param filename: input CSV file
    :param output_dir: output directory where to put output files
    :param grouping_criteria: which group "column" to use as the groups to calculate relatedness by
    :return:  all the group objects
    """
    global analysis_name
    groups, individuals = load_data(filename)
    missing_indivs = link_individuals(individuals)
    while len(missing_indivs) > 0:
        missing_indivs = link_individuals(missing_indivs)
    check_known_mothers(individuals)
    global missing_mother_fill_rate
    if missing_mother_fill_rate:
        fill_in_missing_mothers(individuals)
    do_traces(groups)
    output_stats(groups, grouping_criteria,
                 os.path.join(output_dir, 'summary_intragroup_relatedness_%s.csv' % analysis_name))
    output_stats_individually(
        groups, os.path.join(output_dir, 'summary_individual_to_group_relatedness_%s.csv' % analysis_name))
    output_data(groups, grouping_criteria, output_dir)
    # dump_info_maps(groups, output_dir)

    groups.update()
    return groups


def load_data(filename):
    """
    Loads the data from the file into individuals and groups
    :param filename: Name of file in csv format with the data
    :return:
    """
    import re
    print('loading data...')
    # open the file
    file_handle = open(filename, 'r')
    # read the data in csv format
    rows = csv.DictReader(file_handle)
    groups = {}
    individuals = {}
    # check that all the column names are accounted for
    global field_definitions
    for _, field_label in field_definitions.iteritems():
        if field_label not in rows.fieldnames:
            print("CSV file %s is missing the field: %s" % (filename, field_label))
            sys.exit(0)
    # collect the rows into individuals and groups

    for row in rows:
        #
        if (not filter_variable) or row[filter_variable].strip(' ') == '1':
            new_indiv = Individual(row)
            global excludes
            if new_indiv.id in excludes:
                continue
            new_indiv.primary_group = new_indiv.primary_group.replace("TABD_A", "TABD")
            new_indiv.primary_group = new_indiv.primary_group.replace("TABD_B", "TABD")
            new_indiv.primary_group = new_indiv.primary_group.replace("TABD_C", "TABD")
            if new_indiv.id not in individuals:
                individuals[new_indiv.id] = new_indiv
            else:
                existing_indiv = individuals[new_indiv.id]
                # print(new_indiv.id + 'has a duplicate')

                if existing_indiv.primary_group != new_indiv.primary_group:
                    if existing_indiv.id[:4] != existing_indiv.primary_group[:4] \
                            and new_indiv.id[:4] == new_indiv.primary_group[:4]:
                        # print('%s replaced group %s->%s' %
                        # (existing_indiv.id, existing_indiv.primary_group, new_indiv.primary_group))
                        individuals[new_indiv.id] = new_indiv
                else:
                    if not existing_indiv.father_id and new_indiv.id:
                        # print('....replaced father id %s %s' % (existing_indiv.father_id, new_indiv))
                        existing_indiv.father_id = new_indiv.father_id
                    if not existing_indiv.mother_id and new_indiv.id:
                        # print('....replaced mother id %s %s' % (existing_indiv.mother_id, new_indiv))
                        existing_indiv.mother_id = new_indiv.mother_id
                    if not existing_indiv.birth_year_str and new_indiv.birth_year_str:
                        existing_indiv.birth_year_str = new_indiv.birth_year_str

                    if existing_indiv.father_id and new_indiv.father_id \
                            and existing_indiv.father_id != new_indiv.father_id:
                        print('%s id inconsistent father is %s or %s' %
                              (existing_indiv, existing_indiv.father_id, new_indiv.father_id))
                        j = 1
                        new_id = "%s_%d" % (new_indiv.id, j)
                        while new_id in individuals:
                            j += 1
                            new_id = "%s_%d" % (new_indiv.id, j)
                        new_indiv.id = new_id
                        individuals[new_id] = new_indiv
                        new_indiv.has_conflict = True
                        existing_indiv.has_conflict = True
    global group_restriction
    global group2_restriction
    g1_ignored_count = 0
    g2_ignored_count = 0
    for indiv_id in individuals:
        new_indiv = individuals[indiv_id]
        if not group_restriction or (group_restriction and re.match(group_restriction, new_indiv.primary_group)):
            if not group2_restriction or\
                    (group2_restriction and re.match(group2_restriction, new_indiv.secondary_group)):
                if new_indiv.primary_group not in groups.keys():
                    groups[new_indiv.primary_group] = Group(new_indiv.primary_group)
                groups[new_indiv.primary_group].add_individual(new_indiv)
            else:
                g2_ignored_count += 1
        else:
            g1_ignored_count += 1

    print "g1_ignored_count", g1_ignored_count
    print "g2_ignored_count", g2_ignored_count
    # close the file
    file_handle.close()
    return groups, individuals


def generate_new_indiv(indiv_id, sex=''):
    fields = {
        field_definitions['id']: indiv_id,
        field_definitions['sex_str']: sex,
        field_definitions['is_alive_str']: '0',
        field_definitions['father_id']: "",
        field_definitions['mother_id']: "",
        field_definitions['birth_year_str']: "",
        field_definitions['simple_name']: "",
        field_definitions['primary_group']: "",
        field_definitions['secondary_group']: "",
    }
    newindiv = Individual(fields)
    return newindiv


def link_individuals(individuals):
    """
    Links individual objects by immediate relationships (mother/father to son and vice versa)
    :param individuals: the list of individuals to link
    """
    import re
    print('linking individuals...')
    missing_father_count = 0
    missing_mother_count = 0
    dead_males = 0
    dead_females = 0
    missing_indivs = {}
    for indiv_name, indiv in individuals.iteritems():
        if indiv.is_alive == 0:
            dead_males += 1 if indiv.sex == 'm' else 0
            dead_females += 1 if indiv.sex == 'f' else 0
        # add the father to this individual and add this individual to the father's offspring list
        if indiv.father_id:
            if indiv.father_id in individuals:
                if individuals[indiv.father_id].sex != 'f':
                    indiv.set_father(individuals[indiv.father_id])
                    individuals[indiv.father_id].sex = 'm'
                else:
                    if indiv.father_id+'f' in individuals:
                        indiv.set_father(individuals[indiv.father_id+'f'])
                    else:
                        indiv2 = generate_new_indiv(indiv.father_id+"m", sex='m')
                        indiv.set_father(indiv2)
                        missing_indivs[indiv.father_id+"m"] = indiv2

            else:
                if indiv.father_id in missing_indivs:
                    # print('found newly created father %s' % indiv.father_id)
                    try:
                        indiv.set_father(missing_indivs[indiv.father_id])
                        missing_indivs[indiv.father_id].sex = 'm'
                    except IncorrectSexException:
                        pass
                else:
                    # print("can't find father id %s" % indiv.father_id)
                    missing_father_count += 1
                    if not re.match('^.{4}0$', indiv.father_id): ######################data specific#
                        new_father = generate_new_indiv(indiv.father_id, sex='m')
                        indiv.set_father(new_father)
                        missing_indivs[indiv.father_id] = new_father

        # add the mother to this individual and add this individual to the mother's offspring list
        if indiv.mother_id:
            if indiv.mother_id in individuals.keys():
                try:
                    indiv.set_mother(individuals[indiv.mother_id])
                    individuals[indiv.mother_id].sex = 'f'
                except IncorrectSexException:
                    pass

            else:
                if indiv.mother_id in missing_indivs:
                    try:
                        indiv.set_mother(missing_indivs[indiv.mother_id])
                        missing_indivs[indiv.mother_id].sex = 'f'
                    except IncorrectSexException:
                        pass
                else:
                    # print("can't find mother id %s" % indiv.mother_id)
                    missing_mother_count += 1
                    if not re.match('^.{4}0$', indiv.mother_id):
                        new_mother = generate_new_indiv(indiv.mother_id, sex='f')
                        indiv.set_mother(new_mother)
                        missing_indivs[indiv.mother_id] = new_mother

    print ('missing fathers: %d  missing mothers: %d' % (missing_father_count, missing_mother_count))
    print ('dead males: %d  dead females: %d' % (dead_males, dead_females))
    return missing_indivs


def fill_in_missing_mothers(individuals):
    """
    Add missing mothers to the list of individuals if they don't already exist.  Our input data had rows for all fathers
    but not all mothers
    :param individuals:
    :return: mothers that were added
    """
    global field_definitions
    global missing_mother_fill_rate
    import random

    k = 0
    new_mothers = {}
    for iname, indiv in individuals.iteritems():
        outstr = ""
        if len(indiv.offspring_list) > 1:
            do_not_have_mothers = []
            existing_mother = None
            for indiv2 in indiv.offspring_list:
                if indiv2.mother:
                    existing_mother = indiv2.mother
                else:
                    do_not_have_mothers.append(indiv2)
            if len(do_not_have_mothers) > 0 and random.random() < missing_mother_fill_rate:
                mother = existing_mother
                if not mother:
                    k += 1
                    new_mother_name = "mother%d: " % k
                    mother = generate_new_indiv(new_mother_name, sex='f')
                    new_mothers[mother.id] = mother
                    outstr = ("Created mother: %s" % new_mother_name)
                for indiv2 in do_not_have_mothers:
                    indiv2.set_mother(mother)
                    indiv2.mother_id = mother.id
                    outstr += "%s (%s), " % (indiv2.simple_name, indiv2.id)

                print outstr
    individuals.update(new_mothers)
    print '%d mothers added' % len(new_mothers)


def do_traces(groups):
    """
    Follows paths within a group
    :param groups:
    """
    print('searching relatedness...')
    for lname, group in sorted(groups.iteritems(), key=lambda x: x[1].smart_name):
        print('\tdoing %s' % group)
        group.do_traces()


def output_data(groups, grouping_criteria, output_dir):
    """
    Print out the cross-charts to csv file
    :param groups: the group objects
    :param grouping_criteria: the grouping criteria to use for intragroup relatedness
    :param output_dir:
    """
    for primary_group, group in sorted(groups.iteritems(), key=lambda x: x[1].smart_name):

        # print('\tsaving %s' % primary_group)
        group.output_relatedness_chart(
            os.path.join(output_dir, 'group_%s-%s_X_chart.csv' %
                         (grouping_criteria,  primary_group.replace(":", "_").replace("?", "_"))), 'mf')
        group.output_relatedness_info_chart(
            os.path.join(output_dir, 'debug-group_%s-%s_X_chart.csv' %
                         (grouping_criteria,  primary_group.replace(":", "_").replace("?", "_"))), 'mf')


def output_stats(groups, grouping_criteria, filename):
    """
    Prints out the calculated stats to the output file
    :param groups:
    :param grouping_criteria:
    :param filename:
    :return:
    """
    global analysis_name
    file_handle = open(filename, 'w')
    file_handle.write(
        'Group,GroupID#,Pedigree_Size,#Alive Males, #Alive Females,'
        'R_mf,SD_mf,R_m,SD_m,R_f,SD_f,R_m/R_f,R_m/R_mf,R_f/R_mf,R_father,R_mother\n')
    r_mfs = []
    r_ms = []
    r_fs = []
    r_mothers = []
    r_fathers = []
    for primary_group, group in sorted(groups.iteritems(),  key=lambda x: x[1].smart_name):
        if len(group.id) < 1:
            continue

        try:

            # Calculates the relatedness of a male to his group
            """
                Relatedness is R = sum [ (1/2)^L(p) ] where p enumerates all paths connecting B and C with unique common
                ancestors and L(p) is the length of path p.
            """
            (mf_rel, mf_stdev) = group.calculate_relatedness('mf')
            (m_rel, m_stdev) = group.calculate_relatedness('m')
            (f_rel, f_stdev) = group.calculate_relatedness('f')
            mother_rel, _ = group.calculate_relatedness('m', through_parents=True)
            father_rel, _ = group.calculate_relatedness('f', through_parents=True)
            print("summary::::::::::::::::::::::")
            print("%s:\t%d\t%d\t%d\t%d\t%d\tr_mf=%0.3f\tr_m=%0.3f\tr_f=%0.3f"
                  % (primary_group, group.size(), group.num_total_males, group.num_alive_males(),
                     group.num_total_females, group.num_alive_females, mf_rel, m_rel, f_rel, ))

            # Handling for missing or incalculable values
            if not (math.isnan(mf_rel) or math.isnan(m_rel) or math.isnan(f_rel) or mf_rel < 0.00001):
                r_mfs.append(mf_rel)
                r_ms.append(m_rel)
                r_fs.append(f_rel)
                r_mothers.append(mother_rel)
                r_fathers.append(father_rel)

                if 0 == mf_rel:
                    r_f_d_mf = float('nan')
                    r_m_d_mf = float('nan')
                else:
                    r_f_d_mf = f_rel / mf_rel
                    r_m_d_mf = m_rel / mf_rel
                if 0 == f_rel:
                    r_m_d_f = float('nan')
                else:
                    r_m_d_f = m_rel/f_rel

                file_handle.write("%s,%s,%d,%d,%d,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f, %0.4f, %0.4f,%0.4f,%0.4f\n"
                                  % (group, group.id, group.size(), group.num_alive_males(),
                                     group.num_alive_females, mf_rel, mf_stdev, m_rel, m_stdev, f_rel, f_stdev,
                                     r_m_d_f, r_m_d_mf, r_f_d_mf, father_rel, mother_rel))

        except NotEnoughGroupMembersException:
            print('group not big')
    # Calculating overall means
    mean_r_mfs = mean(r_mfs)
    mean_r_ms = mean(r_ms)
    mean_r_fs = mean(r_fs)
    mean_r_fathers = mean(r_fathers)
    mean_r_mothers = mean(r_mothers)

    file_handle.write("%s,%s,%d,%d,%d,%0.3f,%0.3f, %0.3f,%0.3f,%0.3f,%0.3f,%0.3f, %0.3f, %0.3f, %0.3f, %0.3f\n"
          % ('avg', '', 0, 0, 0,  mean_r_mfs, stdev(r_mfs), mean_r_ms, stdev(r_ms), mean_r_fs, stdev(r_fs),
                      mean_r_ms/mean_r_fs, mean_r_ms/mean_r_mfs, mean_r_fs/mean_r_mfs, mean_r_fathers, mean_r_mothers))
    print 'Relatedness by %s (%s):' % (grouping_criteria, analysis_name)
    print 'mean relatedness through males: %0.3f' % mean_r_ms
    print 'mean relatedness through females: %0.3f' % mean_r_fs
    print 'mean relatedness through both: %0.3f' % mean_r_mfs
    # print 'mean relatedness through mother: %0.3f' % mean_r_mothers
    # print 'mean relatedness through father: %0.3f' % mean_r_fathers
    file_handle.close()


def output_stats_individually(groups, filename):
    """
    Calculates individual stats (for debugging purposes)
    :param groups:
    :param filename:
    :return:
    """
    file_handle = open(filename, 'w')
    file_handle.write(
        'Indiv_ID,Indiv_Name,Group,GroupID,Pedigree_Size,Male_Pop_Size,'
        'R_mf,R_m,R_f,Father,Mother,PGF,PGM,MGF,MGM,isAlive\n')
    for primary_group, group in sorted(groups.iteritems(), key=lambda x: x[1].smart_name):
        if len(group.alive_males) < 5:
            continue
        for indiv in group.alive_males:
            try:
                mf_rel = indiv.calculate_relatedness_to_list(group.alive_males, 'mf')
                m_rel = indiv.calculate_relatedness_to_list(group.alive_males, 'm')
                f_rel = indiv.calculate_relatedness_to_list(group.alive_males, 'f')

                pgf = ""
                pgm = ""
                mgf = ""
                mgm = ""
                motherm = ""
                fatherm = ""
                if indiv.father:
                    pgf = indiv.father.father_id
                    pgm = indiv.father.mother_id
                    fatherm = indiv.father.full_moniker()
                if indiv.mother:
                    mgf = indiv.mother.father_id
                    mgm = indiv.mother.mother_id
                    motherm = indiv.mother.full_moniker()

                file_handle.write("%s,%s,%s,%s,%d,%d,%0.3f,%0.3f,%0.3f,%s,%s,%s,%s,%s,%s,%s\n" % (
                    indiv.id, indiv.simple_name, group.id, primary_group, group.size(), group.num_alive_males(),
                    mf_rel, m_rel, f_rel, fatherm, motherm, pgf, pgm, mgf, mgm,  indiv.is_alive))
            except NotEnoughGroupMembersException:
                print('Group not big enough')

    file_handle.close()


def dump_info_maps(groups, output_dir):
    """
    Printout a map of individuals and their data (for debugging)
    :param groups:
    :param output_dir:
    :return:
    """
    import os.path
    for primary_group, group in sorted(groups.iteritems(), key=lambda x: x[1].smart_name):
        dump_info_map(group, os.path.join(output_dir, "zz_relatednessinfo_%s.csv" % group.name.replace(":", "_")))


def dump_info_map(group, filename):
    file_handle = open(filename, 'w')
    print filename
    file_handle.write('indivID, fatherID, motherID, PGF, PGM, MGF, MGM, member, gender, alive\n')
    for indiv in group.alive_males:
        pgf = ""
        pgm = ""
        mgf = ""
        mgm = ""
        motherm = ""
        fatherm = ""
        if indiv.father:
            pgf = indiv.father.father_id
            pgm = indiv.father.mother_id
            fatherm = indiv.father.full_moniker()
        if indiv.mother:
            mgf = indiv.mother.father_id
            mgm = indiv.mother.mother_id
            motherm = indiv.mother.full_moniker()

        file_handle.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%d\n" % (
            indiv.full_moniker(), fatherm,
            motherm, pgf, pgm, mgf, mgm,  indiv.primary_group, indiv.sex, indiv.is_alive))
    file_handle.close()


def mean(array):
    """
    calculate a mean from a list of numbers
    :param array:
    :return:
    """
    if len(array):
        float_nums = [float(x) for x in array]
        return sum(float_nums) / len(array)
    else:
        return float('nan')


def stdev(number_list):
    """
    Calculate the standard deviation from a list of numbers
    :param number_list:
    :return:
    """
    if not len(number_list):
        return float('nan')
    mean1 = mean(number_list)
    b = []
    float_numbers = [float(x) for x in number_list]
    for num in float_numbers:
        b.append((num - mean1) ** 2)
    return math.sqrt(sum(b) / len(b))


def get_smart_name(name_string):
    """
    Create a sortable name
    DUND2 -> DUND00002
    :param name_string:
    :return:
    """
    import re
    m = re.search("\d", name_string)
    if m:
        return "%s%s" % (name_string[:m.start()], name_string[m.start():].rjust(7, "0"))
    else:
        return name_string


class Individual:
    """
    An individual object with associated parameters
    """
    def __init__(self, fields):
        self.father = None
        self.mother = None
        self.offspring_list = []
        self.id = None
        self.father_id = None
        self.mother_id = None
        self.primary_group = None
        self.is_alive = False
        self.sex = 'u'
        self.sex_str = None
        self.is_alive_str = None
        self.traces = {}
        self.birth_year = -1
        self.birth_year_str = ""
        self.simple_name = ""
        self.has_conflict = False
        self.primary_group = ""
        self.secondary_group = ""

        global field_definitions
        for k, v in field_definitions.iteritems():
            setattr(self, k, fields[v].strip(' '))

        self.is_alive = (self.is_alive_str.strip() == "1")
        if self.sex_str == '1' or self.sex_str.lower() == "m":
            self.sex = 'm'
        elif self.sex_str == '0' or self.sex_str.lower() == "f":
            self.sex = 'f'
        else:
            self.sex = ''

        try:
            self.birth_year = int(self.birth_year_str)

        except ValueError:
            pass

        self.smart_name = get_smart_name(self.id)

    def full_moniker(self):
        return "%s (%s)" % (self.id, self.simple_name)

    def set_father(self, father):
        """
        Links individual to the father Individual object
        :param father:
        :return:
        """
        if father.sex == 'f':
            error_msg = '%s (%s) Can''t be a biological father for %s!' % (father.id, father.sex,  self.id)
            print(error_msg)
            raise(IncorrectSexException(error_msg))

        self.father = father
        father.set_offspring(self)

    def set_mother(self, mother):
        """
        Links individual to their mother Individual object
        :param mother:
        :return:
        """
        if mother.sex == 'm':
            raise(IncorrectSexException('%s Can''t be a biological mother! for %s' % (mother.id, self.id)))
        self.mother = mother
        mother.set_offspring(self)

    def set_offspring(self, offspring):
        """
        Links this individual to his/her offspring
        :param offspring:
        :return:
        """
        self.offspring_list.append(offspring)

    def __str__(self):
        return self.id

    def __eq__(self, other):
        if isinstance(other, Individual):
            return self.id == other.id
        else:
            return NotImplemented

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def __hash__(self):
        return self.id.__hash__()

    def start_trace(self):
        """
        Starts a trace looking for all possible known ancestors in common between an individual and other individuals
        :return:
        """
        self.traces = {}
        self.do_trace(Trace(self), self.traces)

    def do_trace(self, current_trace, traces, went_up=True):
        """
        Continues a trace looking for all possible known ancestors (heads) and then all possible related individuals
        through that ancestor
        :param current_trace:
        :param traces: a dictionary traces[individuals traced to][heads]
        :param went_up:
        :return:
        """

        # have i already reached this person before?
        if self not in current_trace.links:
            # no i have not. continue the trace.
            current_trace.add_link(self, went_up)

            # am i at the head of the trace so far?  (the most distant common ancestor?)
            if went_up:
                # yes i am.
                current_trace.head = self

            # if this trace between this individual and the currently found relation
            # is shorter than a previous trace, and both have the same head than replace
            # that trace with this one  (or if there is no current trace between the two)

            # is there is a previous trace to this person?
            if self not in traces:
                # if so, add this trace to the existing traces to this person
                traces[self] = {}
                traces[self][current_trace.head] = (copy.copy(current_trace))
            else:
                # okay, if there is a trace to the person already, does that trace have the same common ancestor?
                if current_trace.head in traces[self]:
                    # so there is a trace with the same common ancestor. but is this the shorter of the two?
                    if current_trace.get_length() < traces[self][current_trace.head].get_length():
                        # it is shorter, then replace with this trace.
                        traces[self][current_trace.head] = (copy.copy(current_trace))
                    else:
                        # if not shorter, don't do anything.
                        pass
                else:
                    # there is no trace with the same common ancestor already, so just add already
                    traces[self][current_trace.head] = (copy.copy(current_trace))
            # trace up to the parents and continue
            if went_up:
                if self.father:
                    new_trace = copy.copy(current_trace)
                    self.father.do_trace(new_trace, traces)

                if self.mother:
                    new_trace = copy.copy(current_trace)
                    self.mother.do_trace(new_trace, traces)

            # trace down to the offsprings and continue
            for offspring in self.offspring_list:
                new_trace = copy.copy(current_trace)
                offspring.do_trace(new_trace, traces, went_up=False)

    def get_relatedness_to(self, indiv2, sexes="mf", through_parents=False):
        """
        Calculates the relatedness between two individuals through a trace
        :param indiv2:
        :param sexes:
        :return:
        """
        r = 0.0
        if indiv2 in self.traces:
            for head, trace in self.traces[indiv2].iteritems():

                if through_parents:
                    if trace.first_link_sex in sexes:
                        r += (.5 ** trace.distance())
                else:
                    if head.sex in sexes:
                        r += (.5 ** trace.distance())
        return r

    def get_verbose_relatedness_to(self, indiv2, sexes="mf"):
        """
        Relateness in text format for debugging
        :param indiv2:
        :param sexes:
        :return:
        """
        temp_str = ''
        if indiv2 in self.traces:
            for head, trace in self.traces[indiv2].iteritems():
                if head.sex in sexes:
                    # temp_str += ('%s-%s(%d);' % (head.id, head.sex, trace.distance()))
                    temp_str += (trace.link_string + "; ")
        return temp_str.lstrip(' ')

    def calculate_relatedness_to_list(self, indivs, sexes="mf", through_parents=False):
        """
        Calculates the relatedness between one individual and
        :param indivs:
        :param sexes:
        :return:
        """
        r = 0.0
        if len(indivs) <= 1:
                return 0
        for indiv in indivs:
            if indiv != self:
                r += self.get_relatedness_to(indiv, sexes, through_parents)

        return r / float(len(indivs) - 1)


def check_known_mothers(individuals):
    """
    How many paternal sibs are also maternal sibs
    :param individuals:
    :return:
    """
    stats = []
    for indiv_id, indiv in individuals.iteritems():
        if indiv.sex == "m" and len(indiv.offspring_list) > 1:
            group_stats = []
            for indiv1, indiv2 in itertools.combinations(indiv.offspring_list, 2):
                if indiv1.mother_id and indiv2.mother_id:
                    group_stats.append(1 if indiv1.mother_id == indiv2.mother_id else 0)
            stats.append(1 if all(group_stats) else 0)

    # print stats
    print "number of sets of paternal siblings siblings:", len(stats)
    print "percentage of sets with same mother:", sum(stats) / float(len(stats))


class NotEnoughGroupMembersException(Exception):
    pass


class IncorrectSexException(Exception):
    pass


class Group:
    """
    Group object
    """

    def __init__(self,  primary_group):

        self.indivs = {}
        self.alive_males = []
        self.num_total_males = 0
        self.num_total_females = 0
        self.num_alive_females = 0
        self.id = primary_group
        self.smart_name = get_smart_name(self.id)

    # link individual to this group
    def add_individual(self, indiv):

        self.indivs[indiv.id] = indiv
        global age_restriction
        if indiv.sex == 'm' and indiv.birth_year and \
                (not int(age_restriction) or indiv.birth_year > int(age_restriction)):
            self.num_total_males += 1

        if indiv.is_alive and indiv.sex == 'm':
            if not int(age_restriction) or indiv.birth_year > int(age_restriction):
                self.alive_males.append(indiv)
        if indiv.is_alive and indiv.sex == 'f':
                self.num_alive_females += 1
        if indiv.sex == 'f':
            self.num_total_females += 1

    # return the size of the group
    def size(self):
        return len(self.indivs)

    # return the number alive
    def num_alive_males(self):
        return len(self.alive_males)

    def do_traces(self):
        """
        start doing traces for each living male in the group
        """

        for indiv in self.alive_males:
            indiv.start_trace()

    def __str__(self):
        return self.id

    def calculate_relatedness(self, sexes="mf", through_parents=False):
        """

            Relatedness is R = sum [ (1/2)^L(p) ] where p enumerates all paths connecting B and C with unique common
            ancestors and L(p) is the length of path p.

        """
        relatedness = []

        for indiv1, indiv2 in itertools.combinations(self.alive_males, 2):
            relatedness.append(indiv1.get_relatedness_to(indiv2, sexes, through_parents))
        if len(relatedness) == 0:

            return float('nan'), float('nan')
        return mean(relatedness), stdev(relatedness)

    def output_relatedness_info_chart(self, filename, sexes):
        """
            Create cross charts comparing individuals in a group
        """
        f = open(filename, 'w')
        f.write('name/name,')
        alive_males_sorted = sorted(self.alive_males, key=lambda x: x.smart_name)
        for indiv in alive_males_sorted:
            f.write("%s," % indiv.id)
        f.write('avg.\n')
        for indiv in alive_males_sorted:
            line = indiv.full_moniker() + ","
            for indiv2 in alive_males_sorted:
                line += "%s," % (indiv.get_verbose_relatedness_to(indiv2, sexes))
            line += ('%f\n' % indiv.calculate_relatedness_to_list(alive_males_sorted, sexes))
            f.write(line)
        f.write('avg:,')
        for indiv in alive_males_sorted:
            f.write("%f," % indiv.calculate_relatedness_to_list(alive_males_sorted, sexes))

        f.write('%f' % self.calculate_relatedness(sexes)[0])

        f.close()

    def output_relatedness_chart(self, filename, sexes):
        f = open(filename, 'w')
        f.write('name/name,')
        aliveindivssorted = sorted(self.alive_males, key=lambda x: x.smart_name)
        for indiv in aliveindivssorted:
            f.write("%s," % indiv.simple_name)
        f.write('avg.\n')
        for indiv in aliveindivssorted:
            line = indiv.simple_name + ","
            for indiv2 in aliveindivssorted:
                line += "%f," % (indiv.get_relatedness_to(indiv2, sexes))
            line += ('%f\n' % indiv.calculate_relatedness_to_list(aliveindivssorted, sexes))
            f.write(line)
        f.write('avg:,')
        for indiv in aliveindivssorted:
            f.write("%f," % indiv.calculate_relatedness_to_list(aliveindivssorted, sexes))
        f.write('%f' % self.calculate_relatedness(sexes)[0])

        f.close()


class Trace:
    """
    Trace object, used to track a trace from one individuals to all possible blood relatives.
    """
    def __init__(self, start_indiv, links=None, head=None, link_string=""):
        self.start_indiv = start_indiv
        self.head = None
        self.first_link_sex = 'lkjlk'
        self.links = []
        if head:
            self.head = head
        if links:
            self.links = copy.copy(links)
        self.link_string = link_string

    def __copy__(self):
        return Trace(self.start_indiv, links=self.links, link_string=self.link_string, head=self.head)

    def add_link(self, indiv, went_up):
        if len(self.links) == 0:
            self.first_link_sex = indiv.sex
        self.links.append(indiv)
        direction = "^" if went_up else '!'
        self.link_string += ("%s%s(%s); " % (direction, indiv.id, indiv.sex))

    def __str__(self):
        return self.link_string.replace(self.head, "[%s]" % self.head)

    def distance(self):
        return len(self.links) - 1

    def get_length(self):
        return len(self.links)


main()
