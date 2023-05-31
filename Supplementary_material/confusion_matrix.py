#!/home/giacomo/anaconda3/bin/python
from sys import argv
import numpy as np

def confusion_matrix(input_file,threshold):
    f=open(input_file)

    tp,tn,fp,fn = 0,0,0,0

    for line in f:
        eval = float(line.split()[1])
        classK = int(line.split()[2])
        
        if eval <= threshold and classK == 1: #TRUE POSITIVE
            tp+=1
        elif eval > threshold and classK == 0: #TRUE NEGATIVE
            tn+=1
        elif eval <= threshold and classK == 0: #FALSE POSITIVE
            fp+=1
        elif eval > threshold and classK == 1: #FALSE NEGATIVE
            fn+=1
    f.close()
    cm=np.array([[tp,fp],[fn,tn]])
    acc=(tp+tn)/(tp+tn+fn+fp)
    mcc=((tp*tn)-(fp*fn))/np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    prc=tp/(tp+fp)
    rec=tp/(tp+fn)
    F1=(2*prc*rec)/(prc+rec)
    return cm,acc,mcc,F1


if __name__ == '__main__':
    conf_file=argv[1]
    threshold=float(argv[2])
    cm=confusion_matrix(conf_file,threshold)[0]
    acc=confusion_matrix(conf_file,threshold)[1]
    mcc=confusion_matrix(conf_file,threshold)[2]
    F1=confusion_matrix(conf_file,threshold)[3]
    print("TH:%.1E"%(threshold),"\nCM:",cm,"\nACC:",acc,"\nMCC:",mcc,"\nF1:",F1)
    
