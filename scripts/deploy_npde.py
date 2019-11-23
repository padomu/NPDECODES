#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Created by Filippo Leonardi, edited by Liaowang Huang

import os
import json
import subprocess
import tempfile
import zipfile
import shutil


def copy(file, indir, outdir):
    """
    Aggressively copy file to outdir. File can be a directory and can
    already exists.

    :param file: file (or directory) name
    :param indir: the directory which contains file
    :param outdir: the destination to where file copies to
    :return:
    """
    print("Moving '{}' to '{}'".format(indir + file, outdir + file))
    if file == "":
        shutil.copytree(indir, outdir)
        # Recursively copy an entire directory tree rooted at indir
        # to a directory named outdir and return the destination directory.
    elif file[-1] == "/":
        # if the solution needed to copy is a directory, delete the same directory in outdir (if has)
        # then use copytree to copy the directory
        try:
            shutil.rmtree(outdir + file)  # Delete an entire directory tree
        except:
            pass
        shutil.copytree(indir + file, outdir + file)
    else:
        # copy the file
        shutil.copy2(indir + file, outdir + file)
        """
        shutil.copy(src, dst, *, follow_symlinks=True)
        Copies the file scr to the file or directory dst. src and dst should be strings. If dst specifies a directory, 
        the file will be copied into dst using the base filename from src. Returns the path to the newly created file.
        """


def deploy(filename, indir, outdir, handle_unifdef, with_solution=False, output_name=None):
    """
    Parse a file trough "unifdef" and remove SOLUTIONS ifdef code blocks.
    Must have "unifdef" program. Can toggle if SOLUTIONS is true or false.

    :param filename:  e.g. "1dwaveabsorbingbc.cc"
    :param indir:     "../homeworks/1DWaveAbsorbingBC/mastersolution_tagged"
    :param outdir:    "../homeworks/1DWaveAbsorbingBC/templates"
    :param parse_unifdef: deal with #if ... block if true
    :param with_solution: keep block in #if SOLUTION if true
    :param output_name: if not None, the file after parsing is renamed
    :return:
    """
    # if is_cpp(indir + filename) or is_hpp(indir + filename):
    if handle_unifdef:
        cmd = ["unifdef"]
        # if with_internal:
        #     cmd.append("-DINTERNAL=1")
        # else:
        #     cmd.append("-DINTERNAL=0")
        if with_solution:
            cmd.append("-DSOLUTION=1")
        else:
            cmd.append("-DSOLUTION=0")

        cmd.append(indir + filename)

        if output_name is None:
            print("Preparing {} for ".format(filename) + outdir.split('/')[-2])
        else:
            print("Preparing {} for ".format(output_name) + outdir.split('/')[-2])
        print(" ".join(cmd))

        if output_name is None:
            f = open(outdir + filename, "w")
        else:
            f = open(outdir + output_name, "w")
        subprocess.call(cmd, stdout=f)  # run the command described by cmd
        f.close()
    else:
        copy(filename, indir, outdir)


def mkdir(dir):
    """
    Aggressively create a directory if non existing.
    """
    if os.path.exists(dir):
        shutil.rmtree(dir)
    os.makedirs(dir)


def is_cpp(file):
    """
    Returns True if file looks like a C++ file (header of .cpp)
    """
    # return "cpp" == file.split(".")[-1]
    return file.split(".")[-1] in ["cpp", "cc", "c"]


def is_hpp(file):
    """
    Returns True if file looks like a C++ file (header of .cpp)
    """
    return file.split(".")[-1] in ["hpp", "h"]


def is_py(file_path):
    """
    :param file_path: the path of a file
    :return: true if file is a python file(ends with .py)
    """
    return file_path.split(".")[-1] == "py"


def is_cmake(file_path):
    """
    :param file_path: the path of a file
    :return: true if file with suffix .cmake
    """
    return file_path.split(".")[-1] == "cmake"


def generate_test(problem_dir, problem_dir_source):
    """
    Create test folder in problem_dir (if already exists, delete first)
    Create following files from problem_dir_source/test/,
    problem_dir/test/xxx_test_mastersolution (-DSOLUTION=1)
    problem_dir/test/xxx_test_mysolution (-DSOLUTION=0)

    :param problem_dir: "../homeworks/problem_name/"
    :param problem_dir_source: "../developers/problem_name/"
    :return:
    """
    mkdir(problem_dir + "test/")  # "../homeworks/problem_name/test/"
    fileInTest = os.listdir(problem_dir_source + 'test/')
    for file in fileInTest:
        file_path = problem_dir_source + 'test/' + file
        if is_cpp(file_path):
            # e.g. file = '1dwaveabsorbingbc_test_mastersolution.cc'
            test_problem = file.split('_test_mastersolution')[0]  # test_problem = '1dwaveabsorbingbc'
            suffix = file.split(".")[-1]
            deploy(file, problem_dir_source + 'test/', problem_dir + 'test/', True,
                   with_solution=True, output_name=test_problem + "_test_mastersolution." + suffix)
            deploy(file, problem_dir_source + 'test/', problem_dir + 'test/', True,
                   with_solution=False, output_name=test_problem + "_test_mysolution." + suffix)
        else:
            deploy(file, problem_dir_source + 'test/', problem_dir + 'test/', True,
                   with_solution=True)


def generate_templates_and_mysolution(problem_dir, problem_dir_source):
    """
    First create templates and mysolution directory,
    then copy the files in "mastersolution" in problem_dir_source to them
    (remove "#if SOLUTION" block, that is, -DSOLUTION=0).
    Delete and recreate "mastersolution", copy files with -DSOLUTION=1

    :param problem_dir: "../homeworks/problem_name/"
    :param problem_dir_source: "../developers/problem_name/"
    """
    mkdir(problem_dir + "templates/")   # "../homeworks/problem_name/templates"
    mkdir(problem_dir + "mysolution/")  # "../homeworks/problem_name/mysolution"
    mkdir(problem_dir + "mastersolution/")  # "../homeworks/problem_name/mastersolution"
    fileInMastersolution = os.listdir(problem_dir_source + 'mastersolution/')
    for file in fileInMastersolution:
        file_path = problem_dir_source + 'mastersolution/' + file
        if is_cpp(file_path) or is_hpp(file_path) or is_cmake(file_path) or is_py(file_path):
            handle_unifder = True
        else:
            handle_unifder = False
        deploy(file, problem_dir_source + 'mastersolution/', problem_dir + 'templates/',
               handle_unifder, with_solution=False)
        deploy(file, problem_dir_source + 'mastersolution/', problem_dir + 'mysolution/',
               handle_unifder, with_solution=False)
        deploy(file, problem_dir_source + 'mastersolution/', problem_dir + 'mastersolution/',
               handle_unifder, with_solution=True)


def parse_json(filename):
    """
    Parse a JSON description of the Assignment bundles.
    """
    this_dir = os.path.dirname(os.path.realpath(__file__))
    os.chdir(this_dir)  # change the current working directory to this_dir

    f = open(filename, 'r')

    obj = json.load(f)

    homework_dir = obj["homework_dir"]  # "../homeworks/"
    source_dir = obj["source_dir"]  # "../developers/"
    if not homework_dir.endswith("/"):
        homework_dir += "/"
    if not source_dir.endswith("/"):
        source_dir += "/"

    print("Looking for files in '{}'".format(source_dir))

    for problem in obj["Problems"]:
        if not problem.endswith("/"):
            problem += "/"
        # problem e.g. "1DWaveAbsorbingBC/"
        # first copy the whole folder
        copy(problem, source_dir, homework_dir)
        # then generate the templates, mysolution and test
        problem_dir = homework_dir + problem  # "../homeworks/problem_name/"
        problem_dir_source = source_dir + problem  # "../developers/problem_name/"
        generate_templates_and_mysolution(problem_dir, problem_dir_source)
        if os.path.exists(problem_dir_source + "test/"):
            generate_test(problem_dir, problem_dir_source)
        else:
            print("There is no unit test in {}".format(problem_dir_source))


if __name__ == "__main__":
    try:
        subprocess.call(["unifdef", "--help"])
    except FileNotFoundError:
        print("You must install 'unifdef'.")
        exit(-1)

    parse_json("assignment_list.json")


