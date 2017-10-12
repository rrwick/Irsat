#!/usr/bin/env python


# Copyright 2015 Ryan Wick

# This file is part of Irsat.

# Irsat is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# Irsat is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.

# You should have received a copy of the GNU General Public License along with
# Irsat.  If not, see <http:# www.gnu.org/licenses/>.

import sys
import subprocess
import os
import argparse
import shutil
import datetime
import ConfigParser

outDir = ""
args = ""
commands = ""
lastContigsFile = ""

def main():
    startTime = datetime.datetime.now()
    global args
    args = vars(getArguments())
    checkArguments()
    printStartMessage()
    readConfigFile()
    makeOutputDirectory()

    # For a new run, the starting iteration is 1,
    # but it can be higher if the user specified
    # a resume.
    startIter = args['r'] + 1
    endIter = startIter + args['i']

    for i in range(startIter, endIter):
        iterStartTime = datetime.datetime.now()
        iterDir = makeIterationDirectory(i)
        printIterationMessage(i)

        buildBowtieIndex(i, iterDir)

        if args['1'] != None and args['2'] != None:
            mapPairedReads(i, iterDir)
        if args['u'] != None:
            mapUnpairedReads(i, iterDir)

        if i > 1:
            addPreviousReads(i, iterDir)

        assemble(i, iterDir)

        if not args['keep']:
            deleteTemporaryDirectories(iterDir)

        iterEndTime = datetime.datetime.now()
        duration = iterEndTime - iterStartTime
        print '   Time to complete iteration:', convertTimeDeltaToReadableString(duration)

    endTime = datetime.datetime.now()
    duration = endTime - startTime
    printFinishedMessage(duration)








def printStartMessage():
    print '\nIrsat: Iterative Read Subset Assembly Tool'
    print '------------------------------------------'






def getArguments():
    parser = argparse.ArgumentParser(description="Irsat: Iterative Read Subset Assembly Tool", add_help=False)

    required = parser.add_argument_group('Required arguments')
    reads = parser.add_argument_group('Read arguments', 'Read files must be given as paired reads in separate files (-1 and -2), unpaired reads in a single file (-u), or both paired and unapired reads (-1, -2 and -u)')
    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument("-h", "--help",
                          action="help", help="show this help message and exit")

    optional.add_argument("-k", "--keep",
                          help="keep all temporary files",
                          action="store_true")

    required.add_argument("-c", metavar="CONFIG",
                          help="Configuration file which specifies commands",
                          required=True)

    required.add_argument("-t", metavar="TARGET",
                          help="FASTA file containing one or more target sequences",
                          required=True)

    required.add_argument("-o", metavar="OUTDIR",
                          help="the output directory",
                          required=True)

    required.add_argument("-i", metavar="ITERATIONS",
                          type=int,
                          help="how many assembly iterations will be run",
                          required=True)

    optional.add_argument("-r", metavar="RESUME",
                          type=int,
                          help="Resume an existing Irsat run starting with this iteration",
                          default=0)

    reads.add_argument("-1", metavar="FIRST",
                       help="file of first reads in pair")

    reads.add_argument("-2", metavar="SECOND",
                       help="file of second reads in pair")

    reads.add_argument("-u", metavar="UNPAIRED",
                       help="file of unpaired reads")


    return parser.parse_args()


def checkArguments():
    global args

    if args['1'] == None and args['2'] == None and args['u'] == None:
        print 'You must specify files for either paired-end reads, unparied reads or both.'
        exit()
    if (args['1'] == None and args['2'] != None):
        print 'If a second mate file is given, then a first mate file is also required.'
        exit()
    if (args['1'] != None and args['2'] == None):
        print 'If a first mate file is given, then a second mate file is also required.'
        exit()

    if args['1'] != None and not os.path.isfile(args['1']):
        print 'The first mate file could not be found.'
        exit()
    if args['2'] != None and not os.path.isfile(args['2']):
        print 'The second mate file could not be found.'
        exit()
    if args['u'] != None and not os.path.isfile(args['u']):
        print 'The unpaired file could not be found.'
        exit()
    if not os.path.isfile(args['t']):
        print 'The target file could not be found.'
        exit()


def readConfigFile():
    global commands
    global args

    configFileName = args['c']

    if not os.path.isfile(configFileName):
        print "\nERROR: The configuration file could not be found.\n"
        exit()

    config = ConfigParser.ConfigParser()
    config.read(configFileName)

    commands = {}

    if config.has_option('Mapping', 'paired reads'):
        commands['map_paired'] = config.get('Mapping', 'paired reads').strip().split()

    if config.has_option('Mapping', 'unpaired reads'):
        commands['map_unpaired'] = config.get('Mapping', 'unpaired reads').strip().split()

    if config.has_option('Assembly', 'paired reads'):
        commands['assemble_paired'] = config.get('Assembly', 'paired reads').strip().splitlines()
        commands['assemble_paired'][:] = [line.split() for line in commands['assemble_paired']]

    if config.has_option('Assembly', 'unpaired reads'):
        commands['assemble_unpaired'] = config.get('Assembly', 'unpaired reads').strip().splitlines()
        commands['assemble_unpaired'][:] = [line.split() for line in commands['assemble_unpaired']]

    if config.has_option('Assembly', 'both'):
        commands['assemble_both'] = config.get('Assembly', 'both').strip().splitlines()
        commands['assemble_both'][:] = [line.split() for line in commands['assemble_both']]

    if config.has_option('Assembly', 'contigs'):
        commands['assemble_contigs'] = config.get('Assembly', 'contigs').strip()
    if config.has_option('Assembly', 'graph'):
        commands['assemble_graph'] = config.get('Assembly', 'graph').strip()


# Commands exist as lists where each component separated by a space is its
# own item.  There are variables in the commands like DIRECTORY that will
# need to be replaced with real values before the command is run.
def replacePartOfCommand(command, find, replace):
    returnCommand = []
    for part in command:
        if (part == find):
            returnCommand.append(replace)
        else:
            returnCommand.append(part)
    return returnCommand

def makeOutputDirectory():
    fullOutPath = os.getcwd() + '/' + args['o'] + '/'

    global outDir
    outDir = os.path.dirname(fullOutPath)

    if not os.path.exists(outDir):
        os.makedirs(outDir)





def makeIterationDirectory(iteration):
    fullIterationPath = getIterationDirectoryFullPath(iteration)
    iterDir = os.path.dirname(fullIterationPath)

    # If the iteration directory already exists, delete
    # it, along with all of its contents.
    if os.path.exists(iterDir):
        shutil.rmtree(iterDir)

    os.makedirs(iterDir)
    return iterDir



def getIterationDirectoryFullPath(iteration):
    global outdir
    iterationStr = "%03d" % iteration
    return outDir + '/' + iterationStr + '/'




def printIterationMessage(iteration):
    print '\nIteration ' + str(iteration) + ": "






# Builds a Bowtie 2 index using the target sequences and
# the most current assembly (n/a for first iteration)
def buildBowtieIndex(iteration, iterDir):
    global args
    global lastContigsFile

    print '   ' + getDateTimeString() + '  Building Bowtie 2 index...',
    sys.stdout.flush()

    # Make a folder for the index
    indexDir = iterDir + '/1_mapping_index'
    os.makedirs(indexDir)

    inputFiles = args['t']
    if iteration > 1:
        inputFiles += ',' + lastContigsFile
    outputFiles = indexDir + '/bowtie2index'
    bowtie2_buildCommand = ['bowtie2-build', inputFiles, outputFiles]

    bowtie2_build = subprocess.Popen(bowtie2_buildCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = bowtie2_build.communicate()

    # DO SOMETHING WITH out AND err HERE?

    if bowtie2_build.returncode != 0:
        print 'Bowtie 2 index construction failed.'
        exit()

    print 'done'



# Use BWA or Bowtie to find reads that either map to the references
# or have a pair that makes to the references.
def mapPairedReads(iteration, iterDir):
    global args
    global commands

    print '   ' + getDateTimeString() + '  Mapping paired reads...',
    sys.stdout.flush()

    # Make a folder for the files
    pairedDir = iterDir + '/2-paired_read_alignments'
    os.makedirs(pairedDir)

    # Prepare file paths
    index = iterDir + '/1_mapping_index/bowtie2index'
    unfilteredBam1 = pairedDir + '/alignments'
    unfilteredBam2 = unfilteredBam1 + '.bam'
    bothBam = pairedDir + '/both.bam'
    justReadBam = pairedDir + '/just_read.bam'
    justMateBam = pairedDir + '/just_mate.bam'
    mergedBam = pairedDir + '/merged.bam'
    filteredReads1 = iterDir + '/filtered_reads_R1.fastq'
    filteredReads2 = iterDir + '/filtered_reads_R2.fastq'

    # Use Bowtie2 to run the alignment
    bowtie2Command = commands['map_paired'][:]
    bowtie2Command = replacePartOfCommand(bowtie2Command, 'INDEX', index)
    bowtie2Command = replacePartOfCommand(bowtie2Command, 'PAIRED_READS_FILE_1', args['1'])
    bowtie2Command = replacePartOfCommand(bowtie2Command, 'PAIRED_READS_FILE_2', args['2'])

    # Use samtools view to compress the sam file to a bam
    samtools_viewCommand = ['samtools', 'view', '-Shu', '-']

    # Use samtools sort to sort the bam file
    samtools_sortCommand = ['samtools', 'sort', '-n', '-', unfilteredBam1]

    # Run the commands!
    bowtie2 = subprocess.Popen(bowtie2Command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    samtools_view = subprocess.Popen(samtools_viewCommand, stdin=bowtie2.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    samtools_sort = subprocess.Popen(samtools_sortCommand, stdin=samtools_view.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    samtools_sort.communicate()

    # Use samtools view to produce three separate files:
    #  -reads containing neither 4 nor 8 (read and mate mapped)
    #  -reads containing 8 but not 4 (read mapped, mate didn't)
    #  -reads containing 4 but not 8 (mate mapped, read didn't)
    samtools_viewCommandBoth = ['samtools', 'view', '-u', '-F', '12', '-o', bothBam, unfilteredBam2]
    samtools_viewBoth = subprocess.Popen(samtools_viewCommandBoth, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    samtools_viewBoth.communicate()

    samtools_viewCommandJustRead = ['samtools', 'view', '-u', '-f', '8', '-F', '4', '-o', justReadBam, unfilteredBam2]
    samtools_viewJustRead = subprocess.Popen(samtools_viewCommandJustRead, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    samtools_viewJustRead.communicate()

    samtools_viewCommandJustMate = ['samtools', 'view', '-u', '-f', '4', '-F', '8', '-o', justMateBam, unfilteredBam2]
    samtools_viewJustMate = subprocess.Popen(samtools_viewCommandJustMate, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    samtools_viewJustMate.communicate()

    if samtools_viewBoth.returncode != 0 or samtools_viewJustRead.returncode != 0 or samtools_viewJustMate.returncode != 0:
        print 'Samtools view filtering failed'
        exit()

    # Merge the BAMs into one file
    samtools_mergeCommand = ['samtools', 'merge', '-n', mergedBam, bothBam, justReadBam, justMateBam]
    samtools_merge = subprocess.Popen(samtools_mergeCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    samtools_merge.communicate()

    # Use bedtools bamtofastq to convert the BAM file to two fastq files
    bamtofastqCommand = ['bedtools', 'bamtofastq', '-i', mergedBam, '-fq', filteredReads1, '-fq2', filteredReads2]
    bamtofastq = subprocess.Popen(bamtofastqCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    bamtofastq.communicate()

    if bamtofastq.returncode != 0:
        print 'BAM to FASTQ conversion failed'
        exit()

    print 'done'



# Use BWA or Bowtie to find reads that either map to the references
# or have a pair that makes to the references.
def mapUnpairedReads(iteration, iterDir):
    global args
    print '   ' + getDateTimeString() + '  Mapping unpaired reads...',
    sys.stdout.flush()

    # Make a folder for the files
    unpairedDir = iterDir + '/2-unpaired_read_alignments'
    os.makedirs(unpairedDir)

    # Prepare file paths
    index = iterDir + '/1_mapping_index/bowtie2index'
    unfilteredBam1 = unpairedDir + '/alignments'
    unfilteredBam2 = unfilteredBam1 + '.bam'
    filteredBam = unpairedDir + '/filtered.bam'
    filteredReads = iterDir + '/filtered_reads_U.fastq'

    # Use Bowtie2 to run the alignment
    bowtie2Command = commands['map_unpaired'][:]
    bowtie2Command = replacePartOfCommand(bowtie2Command, 'INDEX', index)
    bowtie2Command = replacePartOfCommand(bowtie2Command, 'UNPAIRED_READS_FILE', args['u'])

    # Use samtools view to compress the sam file to a bam
    samtools_viewCommand = ['samtools', 'view', '-Shu', '-']

    # Use samtools sort to sort the bam file
    samtools_sortCommand = ['samtools', 'sort', '-n', '-', unfilteredBam1]

    # Run the commands!
    bowtie2 = subprocess.Popen(bowtie2Command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    samtools_view = subprocess.Popen(samtools_viewCommand, stdin=bowtie2.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    samtools_sort = subprocess.Popen(samtools_sortCommand, stdin=samtools_view.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    samtools_sort.communicate()

    # Use samtools view to filter out reads that didn't align
    samtools_viewCommand = ['samtools', 'view', '-u', '-F', '4', '-o', filteredBam, unfilteredBam2]
    samtools_view = subprocess.Popen(samtools_viewCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    samtools_view.communicate()

    if samtools_view.returncode != 0:
        print 'Samtools view filtering failed'
        exit()

    # Use bedtools bamtofastq to convert the BAM file to two fastq files
    bamtofastqCommand = ['bedtools', 'bamtofastq', '-i', filteredBam, '-fq', filteredReads]
    bamtofastq = subprocess.Popen(bamtofastqCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    bamtofastq.communicate()

    if bamtofastq.returncode != 0:
        print 'BAM to FASTQ conversion failed'
        exit()

    print 'done'






# For all iterations after the first, look at the previous iteration's
# filtered reads.  Any that aren't included in this iteration should be added
# so the read set always grows, never shrinks.
def addPreviousReads(iteration, iterDir):
    readsMate1 = iterDir + '/filtered_reads_R1.fastq'
    readsMate2 = iterDir + '/filtered_reads_R2.fastq'
    readsUnpaired = iterDir + '/filtered_reads_U.fastq'
    
    previousIterDir = getIterationDirectoryFullPath(iteration - 1)
    previousReadsMate1 = previousIterDir + 'filtered_reads_R1.fastq'
    previousReadsMate2 = previousIterDir + 'filtered_reads_R2.fastq'
    previousReadsUnpaired = previousIterDir + 'filtered_reads_U.fastq'
    
    # Paired reads
    if args['1'] != None and args['2'] != None:
        addReadsFromOneFileToAnother(previousReadsMate1, readsMate1)
        addReadsFromOneFileToAnother(previousReadsMate2, readsMate2)

    # Unpaired reads
    if args['u'] != None:
        addReadsFromOneFileToAnother(previousReadsUnpaired, readsUnpaired)
    


# This function adds any reads that are in the first file but not in the
# second file to the second file.  I.e. it merges the read files into the
# second, without any repeats.
def addReadsFromOneFileToAnother(sourceFile, destinationFile):
    readsAlreadyInDestination = makeDictionaryOfReadNames(destinationFile)
    
    source = open(sourceFile, 'r')
    destination = open(destinationFile, 'a')

    addToDestination = False
    for line in source:
        if line[0] == '@':
            if line.strip() not in readsAlreadyInDestination:
                addToDestination = True
            else:
                addToDestination = False
        if addToDestination:
            destination.write(line)




# This function looks at all reads in a FASTQ file (the lines that start with
# @) and stores their names in a dictionary.
def makeDictionaryOfReadNames(fastqFile):
    returnDict = {}
    fastq = open(fastqFile, 'r')

    for line in fastq:
        if line[0] == '@':
            readName = line.strip()
            returnDict[readName] = ""

    return returnDict
    







# Assemble the filtered reads
def assemble(iteration, iterDir):
    global args
    global commands
    global lastContigsFile

    print '   ' + getDateTimeString() + '  Assembling...',
    sys.stdout.flush()

    # Prepare file paths
    readsMate1 = iterDir + '/filtered_reads_R1.fastq'
    readsMate2 = iterDir + '/filtered_reads_R2.fastq'
    readsUnpaired = iterDir + '/filtered_reads_U.fastq'

    # Make a folder for the assembly
    assemblyDir = iterDir + '/3-assembly'
    os.makedirs(assemblyDir)

    # Paired reads only
    if args['1'] != None and args['2'] != None and args['u'] == None:
        assemblyCommand = commands['assemble_paired'][:]

    # Unpaired reads only
    elif args['1'] == None and args['2'] == None and args['u'] != None:
        assemblyCommand = commands['assemble_unpaired'][:]

    # Both paired and unpaired reads
    else:
        assemblyCommand = commands['assemble_both'][:]

    # Replace variables in the assembly command.
    assemblyCommand[:] = [replacePartOfCommand(line, 'DIRECTORY', assemblyDir) for line in assemblyCommand]
    assemblyCommand[:] = [replacePartOfCommand(line, 'PAIRED_READS_FILE_1', readsMate1) for line in assemblyCommand]
    assemblyCommand[:] = [replacePartOfCommand(line, 'PAIRED_READS_FILE_2', readsMate2) for line in assemblyCommand]
    assemblyCommand[:] = [replacePartOfCommand(line, 'UNPAIRED_READS_FILE', readsUnpaired) for line in assemblyCommand]

    # Execute each line of the assembly commands.  For some assemblers, this
    # may only be one line.  Others, like Velvet, may have multiple lines.
    for command in assemblyCommand:
        assemblyProcess = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = assemblyProcess.communicate()

    # Copy the contigs file to the iteration directory
    contigsFile = assemblyDir + '/' + commands['assemble_contigs']
    if os.path.isfile(contigsFile):
        contigsDestination = iterDir + '/' + os.path.basename(contigsFile)
        shutil.copyfile(contigsFile, contigsDestination)
        lastContigsFile = contigsDestination
    else:
        print "\n\nERROR: the following contig file could not be found:"
        print contigsFile
        print "The configuration file may have the wrong filename or location."
        exit()

    # If a graph was specified in the config file, copy that to the iteration
    # directory too.
    if 'assemble_graph' in commands and commands['assemble_graph'] != "":
       graphFile = assemblyDir + '/' + commands['assemble_graph']
       if os.path.isfile(graphFile):
           graphDestination = iterDir + '/' + os.path.basename(graphFile)
           shutil.copyfile(graphFile, graphDestination)
       else:
           print "\n\nERROR: the following graph file could not be found:"
           print graphFile
           print "The configuration file may have the wrong filename or location."
           exit()

    print 'done'



def deleteTemporaryDirectories(iterDir):
    indexDir = iterDir + '/1_mapping_index'
    pairedDir = iterDir + '/2-paired_read_alignments'
    unpairedDir = iterDir + '/2-unpaired_read_alignments'
    assemblyDir = iterDir + '/3-assembly'

    shutil.rmtree(indexDir)
    if os.path.exists(pairedDir):
        shutil.rmtree(pairedDir)
    if os.path.exists(unpairedDir):
        shutil.rmtree(unpairedDir)
    shutil.rmtree(assemblyDir)





def printFinishedMessage(duration):
    print '\nFinished!'
    print 'Total time to complete:', convertTimeDeltaToReadableString(duration)



def getDateTimeString():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def convertTimeDeltaToReadableString(timeDelta):
    seconds = timeDelta.seconds
    hours = timeDelta.days * 24
    hours += seconds // 3600
    seconds = seconds % 3600
    minutes = seconds // 60
    seconds = seconds % 60
    seconds += timeDelta.microseconds / 1000000.0
    secondString = "{:.1f}".format(seconds)

    returnString = ""
    if hours > 0:
        return str(hours) + ' h, ' + str(minutes) + ' min, ' + secondString + ' s'
    if minutes > 0:
        return str(minutes) + ' min, ' + secondString + ' s'
    return secondString + ' s'


# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()
