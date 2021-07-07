import os
import yaml
#import cv2
from yaml.loader import SafeLoader
import numpy as np
from shutil import copyfile
from enum import Enum
import pandas as pd
import sys


def getYAML(iternumber, K):
    input_yml_path =r"C:\workspace\ubuntu\DefSLAM\scripts\stereo0_template.yaml"
    output_dir = r"C:/workspace/ubuntu/cache/Iter"+str(iternumber)
    os.makedirs(output_dir, exist_ok =True)
    output_yml_path = output_dir+r"\temp.yaml"
    #input file
    fin = open(input_yml_path, "rt")
    #output file to write the result to
    fout = open(output_yml_path, "wt")
    data= fin.read()
    data = data.replace('%PATH_TO_SAVE%', "\""+output_dir+"\"")
    for i in range(K.shape[0]):
        print('%PARAMETER_'+str(i)+'%')
        data = data.replace('%PARAMETER_'+str(i)+'%', str(K[i]))
    fout.write(data)
    fin.close()
    fout.close()
    return output_dir, output_yml_path

def rms_per_sequence (base_dir):

  # rms, rms_std -> rms and standard deviation
  # angIso_mean, angIso_std -> for the angIso
  # angSfN_mean, angSfN_std -> for the anfSfN
  # fraction_of_matches -> inliers / possible matches (from Matches.txt)

  user_input = base_dir
  directory = os.listdir(user_input)

  #searchstring = input('ErrorGTs')
  df_names = []
  df_names_AngIso = []
  df_names_AngSfN = []
  for file in directory:
      
      
      if file.startswith("ErrorGTs"):
          if df_names == []:
              if not os.stat(user_input+"/"+file).st_size == 0:
                  df = pd.read_csv(user_input+"/"+file, header=None).transpose()
                  file = file.strip('ErrorGTs') 
                  file = file.strip('.txt') 
                  df_names.append(file)
          else:
              if not os.stat(user_input+"/"+file).st_size == 0:
                  df_aux = pd.read_csv(user_input+"/"+file, header=None).transpose()
                  file = file.strip('ErrorGTs') 
                  file = file.strip('.txt') 
                  df_names.append(file)
                  df = pd.concat([df, df_aux])
                  
      elif file.startswith("ErrorAngIso"):
          if df_names_AngIso == []:
              if not os.stat(user_input+"/"+file).st_size == 0:
                  df_AngIso = pd.read_csv(user_input+"/"+file, header=None).transpose()
                  file = file.strip('ErrorAngIso') 
                  file = file.split("-")[0]
                  df_names_AngIso.append(file)
          else:
              if not os.stat(user_input+"/"+file).st_size == 0:
                  df_aux = pd.read_csv(user_input+"/"+file, header=None).transpose()
                  file = file.strip('ErrorAngIso') 
                  file = file.split("-")[0]
                  df_names_AngIso.append(file)
                  df_AngIso = pd.concat([df_AngIso, df_aux])
                  
      elif file.startswith("ErrorAngSfN"):
          if df_names_AngSfN == []:
              if not os.stat(user_input+"/"+file).st_size == 0:
                  df_AngSfN = pd.read_csv(user_input+"/"+file, header=None).transpose()
                  file = file.strip('ErrorAngSfN') 
                  file = file.split("-")[0]
                  df_names_AngSfN.append(file)
          else:
              if not os.stat(user_input+"/"+file).st_size == 0:
                  df_aux = pd.read_csv(user_input+"/"+file, header=None).transpose()
                  file = file.strip('ErrorAngSfN') 
                  file = file.split("-")[0]
                  df_names_AngSfN.append(file)
                  df_AngSfN = pd.concat([df_AngSfN, df_aux])
                  
      

  df['frame'] = df_names
  df['frame'] = df['frame'].astype('int32')
  df = df.sort_values(by='frame', ascending=True)
  df.iloc[:,0:-1] = df.iloc[:,0:-1]*1000
  n_matches = df.count(axis=1)-1     

  #df_AngIso['frame'] = df_names_AngIso
  #df_AngIso['frame'] = df_AngIso['frame'].astype('int32')
  #df_AngIso = df_AngIso.sort_values(by='frame', ascending=True)

  #df_AngSfN['frame'] = df_names_AngSfN
  #df_AngSfN['frame'] = df_AngSfN['frame'].astype('int32')
  #df_AngSfN = df_AngSfN.sort_values(by='frame', ascending=True)

  df_matches = pd.read_csv(user_input+"/Matches.txt", sep = ' ', header=None, names=['frame', 'inliers', 'outliers', 'possibleMatches'])
  df_matches['frame'] = df_matches['frame'].astype('int32')
  df_matches = df_matches.sort_values(by='frame', ascending=True)
  
  rms = np.nanmean(df.iloc[:,0:-1])
  rms_std = np.nanstd(df.iloc[:,0:-1])

  angIso_mean =0 # np.nanmean(df_AngIso.iloc[:,0:-1])
  angIso_std = 0#np.nanstd(df_AngIso.iloc[:,0:-1])

  angSfN_mean =0# np.nanmean(df_AngSfN.iloc[:,0:-1])
  angSfN_std =0# np.nanstd(df_AngSfN.iloc[:,0:-1])

  fraction_of_matches = df_matches['inliers'].sum().sum()/df_matches['possibleMatches'].sum().sum()

  return rms, rms_std, angIso_mean, angIso_std, angSfN_mean, angSfN_std, fraction_of_matches

def runDefSLAM(iternumber, params):
    execution_path = r"C:\workspace\ubuntu\DefSLAM\Apps\Release\DefSLAMGT.exe"
    orb_voc_path = r"C:\workspace\ubuntu\DefSLAM\Vocabulary\ORBvoc.txt"
    output_dir, input_yml_path = getYAML(iternumber, K=params)
    input_path =r"C:\workspace\ubuntu\MandalaDataset\Mandala1"
    input_path_image = os.path.join(input_path + r"\images")
    input_path_time = os.path.join(input_path + r"\timestamps\timestamps_short.txt")
    exe_str = execution_path
    exe_str +=" " + orb_voc_path
    exe_str +=" " + input_yml_path
    exe_str +=" " + input_path_image
    exe_str +=" " + input_path_image
    exe_str +=" " + input_path_time
    print(exe_str)
    os.system('cmd /c "'+ exe_str)
    rms, rms_std, angIso_mean, angIso_std, angSfN_mean, angSfN_std, fraction_of_matches = rms_per_sequence(output_dir)
    return rms
#runDefSLAM(1, np.array([1.0]))    

class TwiddleState(Enum):
    change = 0
    start = 1

def twiddle(eval, K, dK, path=""):
    max_it = 100
    opt_dim = 0
    sizeParam = K.shape[0]
    twiddle_state = TwiddleState.start
    iter_number = 0
    optimization_value_best = eval(iter_number, K)
    K[opt_dim] = K[opt_dim] + dK[opt_dim]
       
    stop_norm = np.linalg.norm(dK)*0.1
    
    while (np.linalg.norm(dK) > stop_norm and iter_number<max_it ):
        iter_number= iter_number+1
        optimization_value = eval(iter_number, K)
        if(path==""):
            print("Current state vector : ", K.T ,", ")
            print("Current derivative: " , dK.T, ", ")
            print("Current best:" , optimization_value_best, ", current: ", optimization_value)
            print("Iteration best:" , iter_number, ", ")
        else:
            original_stdout = sys.stdout
            with open(path, 'a') as f:
                sys.stdout = f
                print("Current state vector : ", K.T ,", ")
                print("Current derivative: " , dK.T, ", ")
                print("Current best:" , optimization_value_best, ", current: ", optimization_value)
                print("Iteration best:" , iter_number, ", ")
                sys.stdout = original_stdout

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
    def testEvalFunction( iter, K):
        return np.linalg.norm(K)
    K=np.array([1.0,2.0,3.0])
    dK=np.array([0.1,0.2,0.1])
    
    twiddle(testEvalFunction, K, dK)
testTwiddle()

def runTwiddleDefSLAM():
    K = np.array([12000.0, 0.7])
    dK = np.array([500.0, 0.1])
    twiddle(runDefSLAM, K, dK,"C:/workspace/ubuntu/cache/twiddleTestOut.txt")
runTwiddleDefSLAM()


        


