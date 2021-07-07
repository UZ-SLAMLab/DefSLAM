import os
import yaml
#import cv2
from yaml.loader import SafeLoader
import numpy as np
from shutil import copyfile
from enum import Enum

def runDefSLAM():
    execution_path = r"C:\workspace\ubuntu\DefSLAM\Apps\Release\DefSLAMGT.exe"
    orb_voc_path = r"C:\workspace\ubuntu\DefSLAM\Vocabulary\ORBvoc.txt"
    input_yml_path =r"C:\workspace\ubuntu\MandalaDataset\stereo0.yaml"
    input_path =r"C:\workspace\ubuntu\MandalaDataset\Mandala1"
    input_path_image = os.path.join(input_path + r"\images")
    input_path_time = os.path.join(input_path + r"\timestamps\timestamps.txt")
    exe_str = execution_path
    exe_str +=" " + orb_voc_path
    exe_str +=" " + input_yml_path
    exe_str +=" " + input_path_image
    exe_str +=" " + input_path_image
    exe_str +=" " + input_path_time
    print(exe_str)
    os.system('cmd /c "'+ exe_str)
runDefSLAM()
#def loadCahngeSave():
#    input_yml_path =r"C:\workspace\ubuntu\MandalaDataset\stereo0.yaml"
#    cache_dir = r"C:\workspace\ubuntu\cache\temp.yaml"
#    copyfile(input_yml_path, cache_dir)
#    fs = cv2.FileStorage(cache_dir, cv2.FILE_STORAGE_READ)
#    root = fs.root()
#    for name in root.keys():
#        print(name)
#       
#    fs.release()

class TwiddleState(Enum):
    change = 0
    start = 1

def twiddle(eval, K, dK):
    max_it = 100
    opt_dim = 0
    sizeParam = len(K)
    twiddle_state = TwiddleState.start
    K[opt_dim] = K[opt_dim] + dK[opt_dim]
    iter_number = 0
    optimization_value_best = eval(K, iter_number)
    
    while (np.linalg.norm(dK) > 1e-5 and iter_number<max_it ):
        iter_number= iter_number+1
        optimization_value = eval(K, iter_number)
        print("Current result: ", K.T ,", ")
        print("Current steps: " , dK.T, ", ")
        print("Current best:" , optimization_value_best, ", current: ", optimization_value)
        print("Iteration best:" , iter_number, ", ")
        if twiddle_state ==  TwiddleState.start:
            if (optimization_value < optimization_value_best):
                optimization_value_best = optimization_value
                dK[opt_dim] = dK[opt_dim] * 1.5
                opt_dim = (opt_dim+1)  % sizeParam
                K[opt_dim] = K[opt_dim] + dK[opt_dim]
            else:
                K[opt_dim] = K[opt_dim] - 2.0*dK[opt_dim]
                twiddle_state = TwiddleState.change
        elif twiddle_state== TwiddleState.change:
            if optimization_value < optimization_value_best:
                optimization_value_best = optimization_value
                dK[opt_dim] = -dK[opt_dim]
            else:
                K[opt_dim] = K[opt_dim] + dK[opt_dim]
                dK[opt_dim] = dK[opt_dim] / 2.0    
            opt_dim =  (opt_dim + 1)  % sizeParam
            twiddle_state = TwiddleState.start
            K[opt_dim] = K[opt_dim] + dK[opt_dim]
def testTwiddle():        
    def testEvalFunction(K, iter):
        return np.linalg.norm(K)
    K=np.array([1.0,2.0,3.0])
    dK=np.array([0.1,0.2,0.1])
    
    twiddle(testEvalFunction, K, dK)
testTwiddle()

        


