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
def det(a):
    """
    :param a: matrics
    :return: determinant of matrics a
    """
    if len(a) and len(a[0]) is 2:
        return a[0][0] * a[1][1] - a[0][1] * a[1][0]
    sum1 = 0
    for j in range(len(a[0])):
        if j % 2 == 0:
            sign = 1
        else:
            sign = -1
        sum1 += sign * a[0][j] * det(minor(a, 0, j))
    return sum1


def minor(b, row, col):
    """
    :param b: matrics
    :param row: row to remove
    :param col: column to remove
    :return: matrics b without row and col
    """
    if row >= len(b) and col >= len(b):
        return b
    c = makeMatrics(len(b) - 1, len(b) - 1)
    x = 0
    y = 0
    for i in range(len(b)):
        for j in range(len(b[0])):
            if i is not row and j is not col:
                c[x][y] = b[i][j]
                if y is len(c[0]) - 1:
                    x += 1
                    y = 0
                else:
                    y += 1
    return c


# ############# LU ##############
# def findU(a, pivoting=0):
#     """
#     :param a: matrics
#     :param pivoting: indicates- 0 if not using pivoting, else if using pivoting
#     :return: U and inverse L matrices
#     """
#     invL = unitMatrics(makeMatrics(len(a), len(a[0])))
#     for row in range(len(a)):
#         j = row + 1
#         if a[row][row] is not 0:
#             while j < len(a):
#                 if pivoting is not 0:
#                     a = swapRow(a, j, checkPivot(a, j, row))
#                 b = elementalMatrics(a, j, row)
#                 a = multMatrics(b, a)
#                 invL = multMatrics(b, invL)
#             j += 1
#         else:
#             while j < len(a):
#                 if a[j][row] is not 0:
#                     a = swapRow(a, row, j)
#                     break
#                 j += 1
#     return a, invL

def findU(a, pivoting=0):
    """
    :param a: matrics
    :param pivoting: indicates- 0 if not using pivoting, else if using pivoting
    :return: U and inverse L matrices
    """
    invL = unitMatrics(makeMatrics(len(a), len(a[0])))
    for row in range(len(a)):
        j = row + 1
        while j < len(a):
            if pivoting is not 0:   # inverseA
                b = swapRow(a, j, checkPivot(a, j, row))
                a = multMatrics(b, a)
                invL = multMatrics(b, invL)
            elif a[row][row] is 0:  # LU
                k = j + 1
                while k < len(a):
                    if a[k][row] is not 0:
                        b = swapRow(a, row, k)
                        a = multMatrics(b, a)
                        invL = multMatrics(b, invL)
                        break
                    k += 1
            b = elementalMatrics(a, j, row)
            a = multMatrics(b, a)
            invL = multMatrics(b, invL)
            j += 1
    return a, invL


def dispatchU(a, index=0, pivoting=0):  # pivoting is for gaussian instability
    """
    :param a: matrics
    :param index: indicates which option to return from func
    :param pivoting: indicates- 0 if not using pivoting, else if using pivoting
    :return: U, invL or both, with pivoting or without
    """
    U, invL = findU(a, pivoting)
    if index is 0:
        return U
    if index is 1:
        return invL
    return U, invL


def findL(a):
    """
    :param a: matrics
    :return: L matrics
    """
    return inverse(dispatchU(a, 1))


# return row of bigger pivot in col to swap rows
def checkPivot(a, i, j):
    """
    :param a: matrics
    :param i: pivot row index
    :param j: pivot col index
    :return: row index of the biggest element under pivot
    """
    for x in range(i + 1, len(a[0])):
        if abs(a[i][j]) < abs(a[x][j]):
            i = x
    return i


def inverse(a, pivoting=0):
    """
    :param a: matrics
    :param pivoting: indicates- 0 if not using pivoting, else if using pivoting
    :return: inverse matrics a
    """
    if det(a) is 0:
        return
    a, matInverse = dispatchU(a, 2, pivoting)
    a, matInverse = oneOnDiagonal(a, matInverse)
    size = len(a[0]) - 1
    while size > 0:
        for i in range(size):
            b = elementalMatrics(a, i, size)
            a = multMatrics(b, a)
            matInverse = multMatrics(b, matInverse)
        size -= 1
    return matInverse


def oneOnDiagonal(a, matInverse):
    """
    :param a: matrics
    :param matInverse: matrics composed from all the elemental operations on matrics a
    :return: a with 1 on the main diagonal, matInverse updated with the new elemental operations
    """
    b = makeMatrics(len(a), len(a[0]))
    b = unitMatrics(b)
    for i in range(len(a[0])):
        b = unitMatrics(b)
        b[i][i] = 1 / a[i][i]
        a = multMatrics(b, a)
        matInverse = multMatrics(b, matInverse)
    return a, matInverse


# ########## range of error ############
def infNorm(a):
    """
    :param a: matrics
    :return: infinity norm of matrics a
    """
    norm = 0
    for i in range(len(a[0])):
        sumRow = 0
        for j in range(len(a)):
            sumRow += abs(a[i][j])
        norm = max(sumRow, norm)
    return norm


def condA(a, invA):
    """
    :param a: matrics
    :param invA: the inverse of matrics a
    :return: cond of matrics a
    """
    return infNorm(a) * infNorm(invA)


def elementalMatrics(a, i, j):
    """
    :param a: matrics
    :param i: row index
    :param j: column index
    :return: elemental matrics to make a[i][j] = 0
    """
    c = makeMatrics(len(a), len(a[0]))
    c = unitMatrics(c)
    c[i][j] = -1 * (a[i][j] / a[j][j])
    return c


def unitMatrics(c):
    """
    :param c: get matrix
    :return: make her to unit matrix
    """
    for x in range(len(c)):
        c[x][x] = 1  # make c a unit matrix
    return c


def multMatrics(a, b):
    """
    :param a: get matrix
    :param b: get matrix
    :return: new matrix of mul between them
    """
    if len(a[0]) is len(b):
        c = makeMatrics(len(a), len(b[0]))
        for row in range(len(a)):
            for col in range(len(b[0])):
                for x in range(len(a)):
                    c[row][col] += (a[row][x] * b[x][col])
        return c
    return None


def makeMatrics(row, col):
    """
    :param row: get rows of matrix
    :param col: get columns of matrix
    :return: return zero matrix
    """
    c = []
    for i in range(row):
        c += [[0] * col]
    return c


def setMatrics():
    """
    :return:  ask size and numbers for matrics
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
    :param a: original matrics
    :param r1: 1st row
    :param r2: row to swap with
    :return: elemental swap matrics
    """
    c = makeMatrics(len(a), len(a[0]))
    c = unitMatrics(c)
    c[r1] = a[r2]
    c[r2] = a[r1]
    return c


def printMat(a):
    """
    :param a: a matrix to print
    :return: prints in matrix format
    """
    # print('\n'.join(['\t'.join(['{:4}'.format(item) for item in row])
    #                  for row in a]))
    for i in range(len(a)):
        print(a[i], end='\n')
    print("---------------")


def gaussianElimination(a, b):
    """
    :param a: get matrix a
    :param b: get result vector
    :return: print to terminal
    """
    invA = inverse(a, 1)
    cond = condA(a, invA)
    print("inverse A = ")
    printMat(invA)
    print("another presetation = ")
    print(invA)
    print("---------------")
    print("cond = " + str(cond))
    print("---------------")
    print("x = ")
    x = multMatrics(invA, b)
    printMat(x)
    print("another presetation = ")
    print(x)
    print("---------------")


def LUdecomposition(a, b):
    """
    :param a: get matrix a
    :param b: get result vector
    :return: prints LU and solution vector x
     """
    U, invL = dispatchU(a, 2)
    L = inverse(invL)
    print("L = ")
    printMat(L)
    print("another presetation = ")
    print(L)
    print("---------------")
    print("U = ")
    printMat(U)
    print("another presetation = ")
    print(U)
    print("---------------")
    print("x = ")
    x = multMatrics(multMatrics(inverse(U), invL), b)
    printMat(x)
    print("another presetation = ")
    print(L)
    print("---------------")


def driver():
    """
    main function
    :return: prints solution
    """
    a = [[2, 1],
         [4, 3]]

    b = [[1],
         [2]]

    if det(a) is 0:
        print("this matrix is singular")
        return

    if len(a) < 4:
        gaussianElimination(a, b)
    else:
        LUdecomposition(a, b)


driver()

# printMat(inverse([[1, 2, 1], [2, 6, 1], [1, 1, 4]]))  # [4.6,-1.4,-0.8],[-1.4,0.6,0.2],[-0.8,0.2,0.4]
# printMat(inverse(dispatchU(a, 0, 1)))  # [4.6,-1.4,-0.8],[-1.4,0.6,0.2],[-0.8,0.2,0.4]
# printMat(inverse([[2, -3, -5], [1, -1, -1], [-1, 3, 5]]))
# 1.0,  0.0,  1.0
# 2.0, -2.5,  1.5
# -1.0,  1.5, -0.5
# a = [[1,1,1,0],[0,3,1,2],[2,3,1,0],[1,0,2,1]]
# a = [[3,-2,4],[1,0,2],[0,1,0]]
# print(det(a))
# LUdecomposition([[11, 25, 3, 6], [2, 7, 8, 6], [1, 14, 7, 32], [25, 9, 45, 12]], [[1], [2], [1], [2]])
# LUdecomposition([[7,3,-1,2], [3,8,1,-4], [-1,1,4,-1], [2,-4,-1,6]], [[1], [1], [1], [1]])
# gaussianElimination([[2,-3,-5], [1,6,1], [10,3,5]], [[2], [1], [3]])
# LUdecomposition([[5,4,2,6], [4,2,1,0], [0,0,1,1], [2,3,1,5]], [[1], [1], [1], [1]])
# gaussianElimination([[2,-3,-5], [1,6,1], [10,3,5]], [[2], [1], [3]])
# gaussianElimination([[2,-5,3], [0,7,-2], [-1,4,1]], [[1], [2], [1]])
