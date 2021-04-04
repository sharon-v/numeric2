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
    a = [[1, 2], [3, 4]]    # try different sizes!!
    b = [[2, 2], [2, 2]]
    size = len(a)
    # make new matrics c
    c = [[0, 0], [0, 0]]
    for row in range(size):
        for col in range(size):
            for x in range(size):
                c[row][col] += (a[row][x] * b[x][col])
    print(c)


test()
