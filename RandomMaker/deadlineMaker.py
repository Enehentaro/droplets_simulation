import math

def virus_deadline(death_param):
    halfLife = 3960.0   #半減期 1.1 h ( = 3960 sec)
    deadline = - math.log(death_param, 2) * halfLife
    return deadline

import random
random_list = [random.random() for i in range(10000)]

deadline_list = list(map(virus_deadline, random_list))
print(sum(x > 3960.0 for x in deadline_list))

with open('deadline.txt', 'w') as f:
    f.write(str(len(deadline_list)) + "\n")
    for d in deadline_list:
        f.write("%s\n" % str(d))
