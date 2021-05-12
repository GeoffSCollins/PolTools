# STEP4: TSR FINDING ALGORITHM from Rocky which I have modified slightly

from collections import defaultdict



def run_step_four(step_three_file, window_size, chromosome_sizes, output_filename):
    fi = open(step_three_file, 'r')

    DDCHRSIZE = chromosome_sizes

    id_to_process = '1'
    lines_to_process = []
    OFILE = open(output_filename, 'w')
    STEP = window_size


    def TSRF(LIST, STEPSIZE):
        k = 0
        TSR = []

        append_last = True

        while repr(k) != repr(len(LIST)):
            DDL = defaultdict(list)
            for f in LIST:
                DATA = f
                TMP = DATA[0].strip() + ':' + DATA[1].strip() + '-' + DATA[2].strip()
                DDL[TMP] = DATA
            REMOVE = set()
            try:
                MAX = LIST[0]
            except IndexError:
                append_last = False
                break


            for l in LIST[1:]:
                if len(range(max(int(l[1]), int(MAX[1])), min(int(l[2]), int(MAX[2])))) > 0 and int(l[4]) > int(MAX[4]):
                    MAX = l
                elif len(range(max(int(l[1]), int(MAX[1])), min(int(l[2]), int(MAX[2])))) < 1:
                    TSR.append(MAX)
                    MAX = l
            TSR.append(MAX)
            REMOVE = [rn[0].strip() + ':' + '-'.join([repr(nr), repr(nr + STEPSIZE)]).strip() for rn in TSR for nr in
                      range((int(rn[1]) + 1) - STEPSIZE, int(rn[2]))]
            k = len(LIST)
            LIST = [[int(DDL[ls][1]), DDL[ls]] for ls in set(DDL.keys()).difference(REMOVE)]
            LIST.sort(key=lambda x: x[0])
            LIST = [sl[1] for sl in LIST]

        if append_last:
            try:
                TSR.append(MAX)
            except UnboundLocalError:
                pass

        TSR = set(['\t'.join(ts).strip() for ts in TSR])
        TSR = [[int(srrt.strip().split('\t')[1]), srrt.strip()] for srrt in list(TSR)]
        TSR.sort(key=lambda x: x[0])
        TSR = [w[1].strip() for w in TSR]
        return TSR


    while True:
        raw_line = fi.readline()
        # check whether we're at eof
        if not raw_line:
            break
        # get line
        line = raw_line.rstrip().split("\t")
        id = line[0].strip()
        strand = line[-1].strip()
        # same id
        if id == id_to_process.strip():
            lines_to_process.append(line)
        # a new id is coming
        else:
            # do something with your lines
            if len(lines_to_process) > 0:
                for tsrf in TSRF(lines_to_process, STEP):
                    CHECKTSRBUNDRY = tsrf.strip().split('\t')
                    if DDCHRSIZE[CHECKTSRBUNDRY[0].strip()] > 0 and int(CHECKTSRBUNDRY[2]) <= DDCHRSIZE[
                        CHECKTSRBUNDRY[0].strip()]:
                        OFILE.write(tsrf.strip() + '\n')
            # now proceed with next id
            id_to_process = id.strip()
            lines_to_process = [line]

    for tsrl in TSRF(lines_to_process, STEP):
        CHECKTSRBUNDRY = tsrl.strip().split('\t')
        if DDCHRSIZE[CHECKTSRBUNDRY[0].strip()] > 0 and int(CHECKTSRBUNDRY[2]) <= DDCHRSIZE[CHECKTSRBUNDRY[0].strip()]:
            OFILE.write(tsrl.strip() + '\n')
    OFILE.close()
