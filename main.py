""" assignment 2
    team members: Hadar Amsalem, Sharon Vazana
    git link: https://github.com/sharon-v/numeric2.git"""


# gaussian elimination--------------
# non singular if inverse exists
# inverse exists if (det != 0)
# solve with (inverse A) * (vector b)
# LU--------------------------------
# singular if no inverse -> (det == 0)
# solve with (inverse U) * (inverse L) * (vector b)

# (m*n) * (n*p) = m*p

# ############# det ##############


# def det(a):
#     if len(a) and len(a[0]) is 2:
#         return a[0][0] * a[1][1] - a[0][1] * a[1][0]
#     i = 0
#     x = 1
#     for j in range(len(a[0])):
#         return a[i][j]


# def minor(b, row, col):
#     if row >= len(b) and col >= len(b):
#         return b
#     c = makeMatrics(len(b) - 1, len(b) - 1)
#     x = 0
#     y = 0
#     for i in range(len(b)):
#         for j in range(len(b[0])):
#             if i is not row and j is not col:
#                 c[x][y] = b[i][j]
#                 if y is len(c[0]) - 1:
#                     x += 1
#                     y = 0
#                 else:
#                     y += 1
#     return c


def findU(a):
    matL = unitMatrics(makeMatrics(len(a), len(a[0])))
    for row in range(len(a)):
        j = row + 1
        if a[row][row] is not 0:
            while j < len(a):
                if a[j][row] is not 0:
                    b = elementalMatrics(a, j, row)
                    a = multMatrics(b, a)
                    b[j][row] *= -1
                    matL = multMatrics(b, matL)
                j += 1
        else:
            while j < len(a):
                if a[j][row] is not 0:
                    a = swapRow(a, row, j)
                    break
                j += 1
        # what to do if all under pivot are zero
    # return a, matL
    printMat(a)
    print("--------------")
    printMat(matL)


def elementalMatrics(a, i, j):
    c = makeMatrics(len(a), len(a[0]))
    for x in range(len(a)):
        c[x][x] = 1  # make c a unit matrix
    c[i][j] = -1 * (a[i][j] / a[j][j])
    return c


def unitMatrics(c):
    for x in range(len(c)):
        c[x][x] = 1  # make c a unit matrix
    return c


# def multMatrics(a, b):
#     c = makeMatrics(len(a), len(b[0]))
#     for row in range(len(a)):
#         for col in range(len(b[0])):
#             for x in range(len(a)):
#                 c[row][col] += (a[row][x] * b[x][col])
#     return c

def multMatrics(b, a):
    c = makeMatrics(len(b), len(a[0]))
    for row in range(len(b)):
        for col in range(len(a[0])):
            for x in range(len(b)):
                c[row][col] += (b[row][x] * a[x][col])
    return c


def makeMatrics(row, col):
    """
    :param row: amount of rows
    :param col: amount of columns
    :return: zero matrix at the arguments size
    """
    c = []
    for i in range(row):
        c += [[0] * col]
    return c


def setMatrics():
    """
    :return: input matrix from user
    """
    row = int(input('Enter rows >>> '))
    col = int(input('Enter columns >>> '))
    c = makeMatrics(row, col)
    for i in range(row):
        for j in range(col):
            c[i][j] = int(input('Enter number: '))
    return c


def swapRow(a, r1, r2):
    """
    :param a: original matrix
    :param r1: first row to swap
    :param r2: the row to swap with
    :return: the original matrix with the rows swapped
    """
    if r2 < len(a) and r1 < len(a):
        temp = a[r1]
        a[r1] = a[r2]
        a[r2] = temp
    return a


def solveLU(invU, invL, b):
    """
    :param invU: inverse U matrix
    :param invL: inverse L matrix
    :param b: solution vector b
    :return: multiplication of the three
    """
    return multMatrics(invU, multMatrics(invL, b))


def printMat(a):
    """
    :param a: a matrix to print
    :return: prints in matrix format
    """
    for i in range(len(a)):
        j = 0
        while j < len(a[0]):
            print(a[i][j], end=" ")
            j += 1
        print("")


def drive():
    findU([[1, 2, 1], [2, 6, 1], [1, 1, 4]])
    # print(minor([[1, 2, 3], [4, 5, 6], [7, 8, 9]], 2, 2))
    c = setMatrics()
    print(c)
    print(swapRow(c, 1, 0))


drive()
