import PySimpleGUI as sg
import sys
import subprocess

sg.change_look_and_feel('Reddit')      # Add some color to the window

# String finder
def longestSubstringFinder(string1, string2):
    answer = ""
    len1, len2 = len(string1), len(string2)
    for i in range(len1):
        match = ""
        for j in range(len2):
            if (i + j < len1 and string1[i + j] == string2[j]):
                match += string2[j]
            else:
                if (len(match) > len(answer)): answer = match
                match = ""
    return answer

# Please check Demo programs for better examples of launchers
def ExecuteCommandSubprocess(command, *args, timeout=None, wait=False):
    sp = subprocess.Popen([command, *args], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = ''
    for line in sp.stdout:
        line = line.decode(errors='replace' if (sys.version_info) < (3, 5) else 'backslashreplace').rstrip()
        output += line
        print(line)
        win1.Refresh() if win1 else None        # yes, a 1-line if, so shoot me
    retval = sp.wait(timeout)
    return (retval, output)

# ------ Column Definition ------ #
column1 = [[sg.Text('General Parameters', justification='left', font=('Arial', 14))],
            # First entry define outDir
            [sg.Text('1. Set output directory:')],
            [sg.Input(key='outdir'), sg.FolderBrowse(key='browse_1')],
            # Number of threads
            [sg.Text('2. Set the number of threads to use:')],
            [sg.InputText('2', key='n_threads')],
            # Output Prefix
            [sg.Text('3. Set output prefix:')],
            [sg.InputText('out', key='outprefix')],
            # Find YAML file
            [sg.Text('4. Set YAML with additional parameters. (OPTIONAL):')],
            [sg.Text('Yaml File'), sg.Input(key='yaml'), sg.FileBrowse(key='browse_2')],
            # Definition of assembly type
            [sg.Text('5. Choose assembly mode (choose only one):')],
            [sg.Checkbox(' Long reads only assembly', size=(10,1), key='longreadsonly'),  sg.Checkbox(' Short reads only assembly', key='shortreadsonly'),
            sg.Checkbox(' Hybrid (directly)', default=True, key='hybrid_1')],
            [sg.Checkbox(' Hybrid assembly (by polishing long reads contigs)', key='hybrid_2')],
            # Working Directory
            [sg.Text('Set working directory:'), sg.Input(key='working_dir'), sg.FolderBrowse(key='browse_3')]
]

column2 = [[sg.Text('Assembly Parameters', justification='right', font=('Arial', 14))],
           # Definition of assembler
           [sg.Text('Choose assembler (you can choose more than one):')],
           [sg.Checkbox(' Canu', key='use_canu'),  sg.Checkbox(' Flye', key='use_flye'),
           sg.Checkbox(' Unicycler', key='use_unicycler'), sg.Checkbox(' SPAdes', key='use_spades')],
           [sg.Text('Long reads Parameters', justification='right', font=('Arial', 12))],

           # Gather long sequencing reads
           [sg.Text('Set input long reads (Nanopore or Pacbio) in fastq):')],
           [sg.Text('Required for Hybrid or long reads only assembly modes')],
           [sg.Input(key='longreads_fastq'), sg.FileBrowse(key='browse_4')],

           # Check long reads technology
           [sg.Text('Is it nanopore or pacbio?')],
           [sg.Checkbox(' Nanopore', key='input_nanopore'), sg.Checkbox(' Pacbio', key='input_pacbio')],

           # Genome size
           [sg.Text('Expected genome size (e.g. 5.6m):')],
           [sg.Text('Required for Canu and Flye')],
           [sg.Input(key='genomeSize')],

           # Input for polishing
           [sg.Text('Dir containing FAST5 (ONT) for polishing long reads only assemblies (Optional)')],
           [sg.Input(key='fast5_dir'), sg.FolderBrowse(key='browse_5')],
           [sg.Text('Dir containing Pacbio bam and bai files for polishing long reads only assemblies (Optional)')],
           [sg.Input(key='pacbio_bam'), sg.FolderBrowse(key='browse_6')],

           # Gather short sequencing reads
           [sg.Text('Short reads Parameters', justification='right', font=('Arial', 12))],
           [sg.Text('Set input short reads (Illumina) in fastq):')],
           [sg.Text('Required for Hybrid or long reads only assembly modes')],
           [sg.Text('Choose unpaired reads'), sg.Input(key='shortreads_unpaired'), sg.FileBrowse(key='browse_7')],
           [sg.Text('Paired reads must be in the same directory and always with .fastq suffix')],
           [sg.Text('Choose read pairs'), sg.Input(key='shortreads_paired'), sg.FilesBrowse(key='browse_8')]
]

# Design home pattern
layout = [[sg.Text('MpGAP - A generic multiplatform genome assembly pipeline (Set parameters before running pipeline)', justification = 'center')],
          [sg.Button('Download Docker images'), sg.Button('General Parameters'),
          sg.Button('Assembly Parameters'), sg.Button('Save Parameters'), sg.Button('Run Pipeline')],
          [sg.Text('Script output....', size=(40, 1))],
          [sg.Output(size=(100, 20))]
]


# Draw window one
win1 = sg.Window('MpGAP', layout)
 ## General Parameters
win2_active=False
 ## Assembly Parameters
win3_active=False

while True:
    ev1, vals1 = win1.Read(timeout=100)
    if ev1 is None:
        break

    if ev1 == 'General Parameters'  and not win2_active:
        win2_active = True
        win1.Hide()
        layout2 = [[sg.Column(column1)],
    [sg.Button('Done')]]

        win2 = sg.Window('General Parameters', layout2)
        while True:
            ev2, vals2 = win2.Read()
            if ev2 is None or ev2 == 'Done':
                # Close
                win2.Close()
                win2_active = False
                win1.UnHide()
                break

    if ev1 == 'Assembly Parameters'  and not win3_active:
        win3_active = True
        win1.Hide()
        layout3 = [[sg.Column(column2)],
    [sg.Button('Done')]]

        win3 = sg.Window('Assembly Parameters', layout3)
        while True:
            ev3, vals3 = win3.Read()
            if ev3 is None or ev3 == 'Done':
                # Close
                win3.Close()
                win3_active = False
                win1.UnHide()
                break
    if ev1 == 'Download Docker images':
        ExecuteCommandSubprocess('docker pull nextflow/nextflow ; docker pull fmalmeida/compgen:BACANNOT', '')
    if ev1 == 'Save Parameters':

## GENERAL PARAMETERS
        # Get file names based on layout indexes
        outDir = ['params.outDir = "', vals2['outdir'], '"']
        paramsOutDir = "".join(outDir)
        # Get threads
        threads = ['params.threads', vals2['n_threads']]
        paramsThreads = " = ".join(threads)
        # Get Prefix
        prefix = ['params.prefix = "', vals2['outprefix'], '"']
        paramsPrefix = "".join(prefix)
        # Get yaml
        yamlFile = ['params.yaml = "', vals2['yaml'], '"']
        paramsYaml = "".join(yamlFile)
        # Get assembly type
        if vals2['longreadsonly'] == True:
            vals2['assembly_type'] = 'longreads-only'
            vals2['illumina_polish_longreads_contigs'] = 'false'
        elif vals2['shortreadsonly'] == True:
            vals2['assembly_type'] = 'shortreads-only'
            vals2['illumina_polish_longreads_contigs'] = 'false'
        elif vals2['hybrid_1'] == True:
            vals2['assembly_type'] = 'hybrid'
            vals2['illumina_polish_longreads_contigs'] = 'false'
        elif vals2['hybrid_2'] == True:
            vals2['assembly_type'] = 'longreads-only'
            vals2['illumina_polish_longreads_contigs'] = 'true'

        illuminaPolish = ['params.illumina_polish_longreads_contigs', vals2['illumina_polish_longreads_contigs'] ]
        paramsIlluminaPolish = " = ".join(illuminaPolish)
        assemblyType = ['params.assembly_type = "', vals2['assembly_type'], '"']
        paramsAssemblyType = "".join(assemblyType)
        # Get working dir
        workDir = vals2['working_dir']
        configFile = [workDir, 'nextflow.config']
        configFile = "/".join(configFile)

## ASSEMBLY PARAMETERS

        paramsTryCanu = 'params.try_canu = false'
        paramsTryFlye = 'params.try_flye = false'
        paramsTryUnicycler = 'params.try_unicycler = false'
        paramsTrySpades = 'params.try_spades = false'

        # Grab assembler
        if vals3['use_canu'] == True:
            paramsTryCanu = 'params.try_canu = true'
        elif vals3['use_flye'] == True:
            paramsTryFlye = 'params.try_flye = true'
        elif vals3['use_unicycler'] == True:
            paramsTryUnicycler = 'params.try_unicycler = true'
        elif vals3['use_spades'] == True:
            paramsTrySpades = 'params.try_spades = true'

        # Get long reads
        lreads = ['params.longreads = "', vals3['longreads_fastq'], '"']
        paramsLreads = "".join(lreads)

        # lreads types
        if vals3['input_nanopore'] == True:
            vals3['lreads_type'] = 'nanopore'
        elif vals3['input_pacbio'] == True:
            vals3['lreads_type'] = 'pacbio'
        else:
            vals3['lreads_type'] = ''

        lreads_type = ['params.lr_type = "', vals3['lreads_type'], '"']
        paramsLreadsType = "".join(lreads_type)

        # Expected Genome Size
        genomeSize = ['params.genomeSize = "', vals3['genomeSize'], '"']
        paramsGenomeSize = "".join(genomeSize)

        # FAST 5 dir
        fast5Dir = ['params.fast5Path = "', vals3['fast5_dir'], '"']
        paramsFast5Dir = "".join(fast5Dir)

        # Pacbio Bam dir
        if vals3['pacbio_bam'] == '':
            pacbio_pattern = ''
        else:
            pacbio_pattern = [vals3['pacbio_bam'], '*.bam']
            pacbio_pattern = "/".join(pacbio_pattern)

        pacbioBamDir = ['params.pacbio_all_bam_path = "', pacbio_pattern, '"']
        paramsPacbioBamDir = "".join(pacbioBamDir)

        # Short reads
        ## Paired
        if vals3['shortreads_paired'] == '':
            pair_pattern = ''
        else:
            read1 = vals3['shortreads_paired'].split(';')[0]
            read2 = vals3['shortreads_paired'].split(';')[1]
            pair_id = longestSubstringFinder(read1, read2)
            pair_pattern = [pair_id, '{1,2}.fastq']
            pair_pattern = "".join(pair_pattern)

        ShortPaired = ['params.shortreads_paired = "', pair_pattern, '"']
        paramsShortPaired = "".join(ShortPaired)
        ## Single
        if vals3['shortreads_unpaired'] == '':
            vals3['shortreads_unpaired'] = ''
        readSingle = ['params.shortreads_single = "', vals3['shortreads_unpaired'], '"']
        paramsShortSingle = "".join(readSingle)

        f = open(configFile, 'w') ;
        print(paramsOutDir, paramsThreads, paramsPrefix, paramsYaml, paramsAssemblyType,
        paramsLreads, paramsIlluminaPolish, paramsTryCanu, paramsTryFlye, paramsTryUnicycler,
        paramsTrySpades, paramsLreadsType, paramsGenomeSize, paramsFast5Dir, paramsPacbioBamDir,
        paramsShortPaired, paramsShortSingle, sep = "\n", file=f)
        f.close()

        print('Configuration file saved in: ', configFile, sep = '')
        print('\nDone! You can now execute the pipeline.\n\n')
        entries = ['echo "nextflow run fmalmeida/MpGAP -c ', configFile, '"']
        cmd = "".join(entries)

    if ev1 == 'Run Pipeline':
        ExecuteCommandSubprocess(cmd, '')
