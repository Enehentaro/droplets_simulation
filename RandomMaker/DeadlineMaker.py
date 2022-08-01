import math
import random

class DeadlineMaker:
    def __init__(self, halfLife:float) -> None:
        self.halfLife = halfLife

    def virus_deadline(self, death_param:float):
        '''
        [0:1]の実数を受け取り、その値に対応する寿命を返す。
        返す値の頻度は半減期に依存。
        '''
        deadline = - math.log(death_param, 2) * self.halfLife
        return deadline

    def get_deadline_list(self, droplets:int) -> list[float]:
        '''
        指定数だけ寿命のリストを返す。
        '''
        random_list = [random.random() for i in range(droplets)]
        deadline_list = list(map(self.virus_deadline, random_list))
        return deadline_list


if __name__ == '__main__':

    halfLife = 3960.0   #半減期 1.1 h ( = 3960 sec)

    dlMaker = DeadlineMaker(halfLife)

    deadline_list = dlMaker.get_deadline_list(10000)
    print(sum(x > 3960.0 for x in deadline_list))

    with open('deadline.txt', 'w') as file:
        file.write(f"{len(deadline_list)}\n")
        for deadline in deadline_list:
            file.write(f"{deadline}\n")
