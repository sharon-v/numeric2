""" assignment 2
    team members: Hadar Amsalem, Sharon Vazana
    git link: https://github.com/sharon-v/numeric2.git"""


# gaussian elimination--------------
# non singular if inverse exists
# inverse exists if (det != 0)
# LU--------------------------------
# singular if no inverse -> (det == 0)
# solve with (inverse U) * (inverse L) * (vector b)

# (m*n) * (n*p) = m*p

# ############# det ##############
def det(A):
    pass


# ########## swap lines ##########


# mult matrics
# def multMatrics(a, b):
#     # assume same size
#     for row in a:
#         for col in a:
#

def test():
    a = [[3, 2], [8, 6]]    # try different sizes!!
    b = [[3, 4], [2, 9]]
    c = makeMatrics(len(a), len(b[0]))
    for row in range(len(a)):
        for col in range(len(b[0])):
            for x in range(len(a)):
                c[row][col] += (a[row][x] * b[x][col])
    print(c)


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


test()
print(setMatrics())
