import glob
import os
import shutil

def copy_csv():

    csv_list = glob.glob("desks_info/*.csv")

    for csv_fname in csv_list:
        bname = os.path.basename(csv_fname)
        office_name = bname.split(".")[0]

        dir_name = f"/run/user/1000/gvfs/smb-share:server=nas1.local,share=enehen_master/konishi/output/droplets/{office_name}_data_droplets"
        if os.path.exists(dir_name):
            shutil.copyfile(
                csv_fname,
                dir_name+"/desks_info.csv"
                )

def get_office_list(path2officeDir:str):
    """

    """
    office_list = []
    for item in os.listdir(path2officeDir):
        if os.path.isdir(os.path.join(path2officeDir, item)):
            office_name = item.split("_")[0]
            office_list.append(office_name)

    return office_list

def get_case_list(path2caseDir:str):
    """

    """
    keys = ["A", "B", "C", "D", "E"]
    case_list = []
    for item in os.listdir(path2caseDir):
        if item.split("_")[-1] in keys:
            case_list.append(item)

    return case_list

if __name__ == '__main__':
    path2officeDir = '/run/user/1000/gvfs/smb-share:server=nas2.local,share=enehen-master2/konishi/output/droplets'
    # office_list = get_office_list(path2officeDir)
    # with open(path2officeDir+'/office_list.txt', mode='w') as f:
    #     f.write('\n'.join(office_list))

    for item in os.listdir(path2officeDir):
        path2caseDir = os.path.join(path2officeDir, item)
        if os.path.isdir(path2caseDir):
            case_list = get_case_list(path2caseDir)
            print(case_list)
            with open(path2caseDir+'/case_list.txt', mode='w') as f:
                f.write('\n'.join(case_list))
            
