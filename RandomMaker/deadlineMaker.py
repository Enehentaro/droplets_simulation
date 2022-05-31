import math
import random

def virus_deadline(death_param):
    halfLife = 3960.0   #半減期 1.1 h ( = 3960 sec)
    deadline = - math.log(death_param, 2) * halfLife
    return deadline

if __name__ == '__main__':

    random_list = [random.random() for i in range(10000)]

    deadline_list = list(map(virus_deadline, random_list))
    print(sum(x > 3960.0 for x in deadline_list))

    with open('deadline.txt', 'w') as file:
        file.write(f"{len(deadline_list)}\n")
        for d in deadline_list:
            file.write(f"{d}\n")
