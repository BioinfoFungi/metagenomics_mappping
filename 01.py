import subprocess as subp
import sys
import bz2
import os

def mybytes(val):
    return bytes(val, encoding='utf-8')
def read_and_split_line(line):
    return line.decode('utf-8').strip().split('\t')
def mapq_filter(marker_name, mapq_value, min_mapq_val):
    if 'GeneID:' in marker_name:
        return True
    else:
        if mapq_value > min_mapq_val:
            return True
    return False

def run_bowtie2(fna_in, outfmt6_out, bowtie2_db, preset, nproc, min_mapq_val, file_format="fasta",
                exe=None, samout=None, min_alignment_len=None, read_min_len=0):
    # checking read_fastx.py
    read_fastx = "read_fastx.py"

    # try:
    #     subp.check_call([read_fastx, "-h"], stdout=DEVNULL, stderr=DEVNULL)
    # except Exception as e:
    #     try:
    #         read_fastx = os.path.join(os.path.join(os.path.dirname(__file__), "utils"), read_fastx)
    #         subp.check_call([read_fastx, "-h"], stdout=DEVNULL, stderr=DEVNULL)
    #     except Exception as e:
    #         sys.stderr.write("OSError: fatal error running '{}'. Is it in the system path?\n".format(read_fastx))
    #         sys.exit(1)
    #
    # # checking bowtie2
    # try:
    #     subp.check_call([exe if exe else 'bowtie2', "-h"], stdout=DEVNULL)
    # except Exception as e:
    #     sys.stderr.write('OSError: "{}"\nFatal error running BowTie2. Is BowTie2 in the system path?\n'.format(e))
    #     sys.exit(1)

    try:
        read_fastx = os.path.join(os.path.join(os.path.dirname(__file__)), read_fastx)
        if fna_in:
            print(" ".join([read_fastx, '-l', str(read_min_len), fna_in]))
            readin = subp.Popen([read_fastx, '-l', str(read_min_len), fna_in], stdout=subp.PIPE, stderr=subp.PIPE)

        else:
            readin = subp.Popen([read_fastx, '-l', str(read_min_len)], stdin=sys.stdin, stdout=subp.PIPE, stderr=subp.PIPE)

        bowtie2_cmd = [exe if exe else 'bowtie2', "--seed", "1992", "--quiet", "--no-unal", "--{}".format(preset),
                       "-S", "-", "-x", bowtie2_db]

        if int(nproc) > 1:
            bowtie2_cmd += ["-p", str(nproc)]

        bowtie2_cmd += ["-U", "-"]  # if not stat.S_ISFIFO(os.stat(fna_in).st_mode) else []

        if file_format == "fasta":
            bowtie2_cmd += ["-f"]
        print(" ".join(bowtie2_cmd))
        p = subp.Popen(bowtie2_cmd, stdout=subp.PIPE, stdin=readin.stdout)
        readin.stdout.close()
        lmybytes, outf = (mybytes, bz2.BZ2File(outfmt6_out, "w")) if outfmt6_out.endswith(".bz2") else (str, open(outfmt6_out, "w"))
        try:
            if samout:
                if samout[-4:] == '.bz2':
                    sam_file = bz2.BZ2File(samout, 'w')
                else:
                    sam_file = open(samout, 'wb')
        except IOError as e:
            sys.stderr.write('IOError: "{}"\nUnable to open sam output file.\n'.format(e))
            sys.exit(1)
        for line in p.stdout:
            if samout:
                sam_file.write(line)

            o = read_and_split_line(line)
            if not o[0].startswith('@'):
                if not o[2].endswith('*'):
                    if (hex(int(o[1]) & 0x100) == '0x0'): #no secondary
                        if mapq_filter(o[2], int(o[4]), min_mapq_val) :  # filter low mapq reads
                            if ((min_alignment_len is None) or
                                    (max([int(x.strip('M')) for x in re.findall(r'(\d*M)', o[5]) if x]) >= min_alignment_len)):
                                outf.write(lmybytes("\t".join([ o[0], o[2].split('/')[0] ]) + "\n"))

        if samout:
            sam_file.close()

        p.communicate()
        read_fastx_stderr = readin.stderr.readlines()
        nreads = None
        avg_read_length = None
        try:
            nreads, avg_read_length = list(map(float, read_fastx_stderr[0].decode().split()))
            if not nreads:
                sys.stderr.write('Fatal error running MetaPhlAn. Total metagenome size was not estimated.\nPlease check your input files.\n')
                sys.exit(1)
            if not avg_read_length:
                sys.stderr.write('Fatal error running MetaPhlAn. The average read length was not estimated.\nPlease check your input files.\n')
                sys.exit(1)
            outf.write(lmybytes('#nreads\t{}\n'.format(int(nreads))))
            outf.write(lmybytes('#avg_read_length\t{}'.format(avg_read_length)))
            outf.close()
        except ValueError:
            sys.stderr.write(b''.join(read_fastx_stderr).decode())
            outf.close()
            os.unlink(outfmt6_out)
            sys.exit(1)

    except OSError as e:
        sys.stderr.write('OSError: "{}"\nFatal error running BowTie2.\n'.format(e))
        sys.exit(1)
    # except IOError as e:
    #     sys.stderr.write('IOError: "{}"\nFatal error running BowTie2.\n'.format(e))
    #     sys.exit(1)
    #
    # if p.returncode == 13:
    #     sys.stderr.write("Permission Denied Error: fatal error running BowTie2."
    #                      "Is the BowTie2 file in the path with execution and read permissions?\n")
    #     sys.exit(1)
    # elif p.returncode != 0:

run_bowtie2('./input/SRS014476-Supragingival_plaque.fasta',
            './output/SRS014476-Supragingival_plaque.fasta.bowtie2out.txt',
            './index/mpa_vJan21_TOY_CHOCOPhlAnSGB_202103.fna',
            'very-sensitive',4,5)