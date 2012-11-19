#!/usr/bin/env python3

if __name__ == '__main__':
    files = []
    for i in range(4):
        files += [open("monitor.csv.%d" % i, 'r')]

    newfile = open("monitor.csv", 'w')

    l = 0
    for line in files[0]:
        if not l:
            newfile.write('# ')
        t0 = line.split()
        newfile.write('%s ' % ' '.join(t0[2:]))
        best = []
        nbindi = []
        avg = []
        if l:
            best += [float(t0[6])]
            nbindi += [float(t0[3])]
            avg += [float(t0[4]) * float(t0[3])]
        for i in range(1,4):
            ti = files[i].readline().split()
            newfile.write('%s ' % ' '.join(ti[3:]))
            if l:
                best += [float(ti[6])]
                nbindi += [float(ti[3])]
                avg += [float(ti[4]) * float(ti[3])]
        if l:
            newfile.write('%d %d' % (max(best), sum(avg) / sum(nbindi)))
        newfile.write('\n')
        l += 1

    print("Done")

