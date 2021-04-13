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
    if len(a) and len(a[0]) is 2:
        return a[0][0] * a[1][1] - a[0][1] * a[1][0]
    sum1 = 0
    for j in range(len(a[0])):
        if j % 2 == 0:
            sign = 1
        else:
            sign = -1
        sum1 += sign*a[0][j]*det(minor(a, 0, j))
    return sum1


def minor(b, row, col):
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


def findU(a):
    matL = unitMatrics(makeMatrics(len(a), len(a[0])))
    a, matL = dispatchU(a, 2)
    for row in range(len(a)):
        j = row + 1
        if a[row][row] is not 0:
            while j < len(a):
                if a[j][row] is not 0:
                    b = elementalMatrics(a, j, row)
                    a = multMatrics(b, a)
                    matL = multMatrics(b, matL)
                j += 1
        else:
            while j < len(a):
                if a[j][row] is not 0:
                    a = swapRow(a, row, j)
                    break
                j += 1
    return a, matL


def dispatchU(a, index=0):
    U, invL = findU(a)
    if index is 0:
        return U
    if index is 1:
        return invL
    return U, invL



# def findLU(a):
#     matL = unitMatrics(makeMatrics(len(a), len(a[0])))
#     for row in range(len(a)):
#         j = row + 1
#         if a[row][row] is not 0:
#             while j < len(a):
#                 if a[j][row] is not 0:
#                     b = elementalMatrics(a, j, row)
#                     a = multMatrics(b, a)
#                     b[j][row] *= -1
#                     matL = multMatrics(b, matL)
#                 j += 1
#         else:
#             while j < len(a):
#                 if a[j][row] is not 0:
#                     a = swapRow(a, row, j)
#                     break
#                 j += 1
#         # what to do if all under pivot are zero
#     return a, matL
#
#
# def dispatchLU(a, index=2):
#     U,L=findLU(a)
#     if index is 0:
#         return U
#     if index is 1:
#         return L
#     else:
#         print(U)
#         print(L)



def findL(a):
    return inverse(findU(a, 1))


def inverse(a):
    if det(a) == 0:
        return
    matInverse = dispatchU(a, 1)
    size = len(a[0])-1
    while size > 0:
        for i in range(size-1):
            b = elementalMatrics(a, i, size)
            a = multMatrics(b, a)
            matInverse = multMatrics(b, matInverse)
        size -= 1
    print(matInverse)
    return oneOnDiagonal(a, matInverse)


def oneOnDiagonal(a, matInverse):
    b = makeMatrics(len(a), len(a[0]))
    for x in range(len(a[0])):
        b[x][x] = 1  # make b a unit matrix
    for i in range(len(a[0])):
        b[i][i] = 1 / a[i][i]
        matInverse = multMatrics(b, matInverse)
    return b


def norma(a):
    sum1 = 0
    norm = 0
    for i in range(len(a[0])):
        for j in range(len(a)):
            sum += abs(a[i][j])
        norm = max(norm, sum1)
    return norm


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


def multMatrics(a, b):
    c = makeMatrics(len(a), len(b[0]))
    for row in range(len(a)):
        for col in range(len(b[0])):
            for x in range(len(a)):
                c[row][col] += (a[row][x] * b[x][col])
    return c


def makeMatrics(row, col):
    c = []
    for i in range(row):
        c += [[0] * col]
    return c


def setMatrics():
    row = int(input('Enter rows >>> '))
    col = int(input('Enter columns >>> '))
    c = makeMatrics(row, col)
    for i in range(row):
        for j in range(col):
            c[i][j] = int(input('Enter number: '))
    return c


def swapRow(a, r1, r2):
    if r2 < len(a) and r1 < len(a):
        temp = a[r1]
        a[r1] = a[r2]
        a[r2] = temp
    return a


def solveLU(invU, invL, b):
    return multMatrics(multMatrics(invU, invL), b)


def drive():
    a = [[1, 2, 1], [2, 6, 1], [1, 1, 4]]
    b = [7, 8, 9]
    if len(a) < 4:
        x = multMatrics(inverse(a), b)
        print("norma: " + norma(x))
    else:
        dispatchU(a)
    # findU([[1, 2, 1], [2, 6, 1], [1, 1, 4]])
    # # print(minor([[1, 2, 3], [4, 5, 6], [7, 8, 9]], 2, 2))
    # c = setMatrics()
    # print(c)
    # print(swapRow(c, 1, 0))
    # print(det([[1, 4, 6], [2, 3, 2], [5, 5, 4]]))


print(inverse(dispatchU(([[1, 2, 1], [2, 6, 1], [1, 1, 4]]), 1)))
# print(inverse([[1, 2, 1], [2, 6, 1], [1, 1, 4]]))
# print(inverse([[1, 4, 6], [2, 3, 2], [5, 5, 4]]))
# drive()
