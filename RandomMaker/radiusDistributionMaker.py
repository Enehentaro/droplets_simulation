import random

def get_radius_list(name:str):
    radius_list = []

    if name == "coughing":
        radius_list.extend([0.45]*386)
        radius_list.extend([0.75]*248)
        radius_list.extend([1.5]*883)
        radius_list.extend([2.5]*1159)
        radius_list.extend([3.5]*1987)
        radius_list.extend([4.5]*2594)
        radius_list.extend([6.25]*36)
        radius_list.extend([8.75]*116)
        radius_list.extend([11.25]*174)
        radius_list.extend([13.75]*528)
        radius_list.extend([31.25]*544)
        radius_list.extend([43.75]*313)
        radius_list.extend([56.25]*212)
        radius_list.extend([68.75]*185)
        radius_list.extend([87.5]*232)
        radius_list.extend([112.5]*151)
        radius_list.extend([187.5]*195)
        radius_list.extend([365]*48)
        radius_list.extend([750]*9)

    elif name == "conversation":
        radius_list.extend([0.45]*363)
        radius_list.extend([0.75]*151)
        radius_list.extend([1.5]*423)
        radius_list.extend([2.5]*1209)
        radius_list.extend([3.5]*1814)
        radius_list.extend([4.5]*2419)
        radius_list.extend([7.5]*907)
        radius_list.extend([12.5]*756)
        radius_list.extend([17.5]*453)
        radius_list.extend([22.5]*302)
        radius_list.extend([35]*726)
        radius_list.extend([47.5]*181)
        radius_list.extend([75]*151)
        radius_list.extend([125]*30)
        radius_list.extend([250]*97)
        radius_list.extend([425]*18)

    return radius_list


if __name__ == '__main__':

    radius_type = "coughing"
    # radius_type = "conversation"

    radius_list = get_radius_list(radius_type)

    print(len(radius_list))

    random.shuffle(radius_list)

    with open("initialRadius_" + radius_type + ".txt", "w") as file:
        file.write(f"{len(radius_list)}\n")
        for r in radius_list:
            file.write(f"{r}\n")