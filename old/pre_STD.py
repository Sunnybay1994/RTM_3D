#!/usr/bin/env python
from numpy import *
import matplotlib.pyplot as plt


def read_src1(fsrc):
    src = []
    data = []
    fp = open(fsrc)
    fline = fp.readline()
    nsrc, nt = fline.split()
    nsrc = int(nsrc)
    nt = int(nt)
    for i in range(nsrc):
        fline = fp.readline()
        ix, iy, iz, Ey = fline.split()
        ix = int(ix)
        iy = int(iy)
        iz = int(iz)
        src.append([ix, iy, iz])
    for i in range(nsrc):
        fline = fp.readline()
        # print fline.split()
        dum = [double(j) for j in fline.split()]
        # print dum
        data.append(dum)
    return nsrc, nt, src, data


def read_src2(fsrc):
    src = []
    data = []
    fp = open(fsrc)
    fline = fp.readline()
    nsrc, nt = fline.split()
    nsrc = int(nsrc)
    nt = int(nt)
    for i in range(nsrc):
        fline = fp.readline()
        ix, iy, iz, Ey = fline.split()
        ix = int(ix)
        iy = int(iy)
        iz = int(iz)
        src.append([ix, iy, iz])
    for i in range(nsrc):
        dum = []
        for j in range(nt):
            fline = fp.readline()
            # print fline.split()
            dum.append(double(fline))
            # print dum
        data.append(dum)
    return nsrc, nt, src, data


def get_src(nsrc, dir1='./Input/', dir2='./RTM/Input/'):
    std_src = []
    std_data = []
    for i in range(nsrc):
        f1 = dir1 + 'src.in_' + str(i + 1).zfill(4)
        f2 = dir2 + 'src.in_' + str(i + 1).zfill(4)
        nsrc1, nt1, src1, data1 = read_src1(f1)
        nsrc2, nt2, src2, data2 = read_src2(f2)
        # print len(data2[0])
        if nsrc1 == 1:
            for j in range(nsrc2):

                if src2[j] == src1[0]:
                    print(j) 
                    std_src.append(src2[j])
                    std_data.append(data2[j])
        else:
            print("nsrc1 != 1") 
    return std_src, std_data


def write_std(std_src, std_data):
    fp = open('src.in_std', 'w')
    nsrc = len(std_data[:, 0])
    nt = len(std_data[0, :])
    fp.write("%d %d\n" % (nsrc, nt))
    for i in range(nsrc):
        fp.write("%d %d %d Ey\n" %
                 (std_src[i][0], std_src[i][1], std_src[i][2]))
    savetxt(fp, std_data)


if __name__ == "__main__":
    std_src, std_data = get_src(289)

    # print len(std_data)
    std_data = array(std_data)
    write_std(std_src, std_data)

    print(shape(std_data)) 
    plt.figure()
    plt.imshow(std_data, cmap='gray', extent=(0, 100, 200, 0))
    plt.savefig('gather_std.png')
