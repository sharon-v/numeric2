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
    for i in range(len(a)):
        if a[i][i] is 0:
            while j < len(a):
                if a[i][j] is not 0:




def multMatrics(a, b):
    c = makeMatrics(len(a), len(b[0]))
    for row in range(len(a)):
        for col in range(len(b[0])):
            for x in range(len(a)):
                c[row][col] += (a[row][x] * b[x][col])


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
    return multMatrics(multMatrics(invU,invL),b)


def drive():
    # print(minor([[1, 2, 3], [4, 5, 6], [7, 8, 9]], 2, 2))
    c = setMatrics()
    print(c)
    print(swapRow(c, 1, 0))


drive()
