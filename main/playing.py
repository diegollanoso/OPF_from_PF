import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# Load the Excel file
#file_path = r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\generator_load_data.xlsx"
#df = pd.read_excel(file_path)

# Display the first few rows of the dataframe
#print(df.head())



#freq = df['Freq'].tolist()
#time = df['Time'].tolist()

#genAGC = df['SumAGC']
#genCPF = df['SumNoAGC']

# Plot genAGC vs time
# plt.figure()
# plt.plot(time, genAGC, label='genAGC')
# plt.xlabel('Time')
# plt.ylabel('genAGC')
# plt.title('genAGC over Time')
# plt.legend()
# plt.savefig('genAGC_over_time.png')
# plt.close()

# # Plot genCPF vs time
# plt.figure()
# plt.plot(time, genCPF, label='genCPF')
# plt.xlabel('Time')
# plt.ylabel('genCPF')
# plt.title('genCPF over Time')
# plt.legend()
# plt.savefig('genCPF_over_time.png')
# plt.close()


# P_recp_list = []

# tstop_cpf = 30
# t_final = 200
# t_step = 5
# n_step = len(range(tstop_cpf+t_step, t_final, t_step))
# tiempo = np.zeros(n_step)
# j=0

# for t_int in range(tstop_cpf+t_step, t_final, t_step):
#     print('t_int = ' + str(t_int))
#     tiempo[j] = t_int
#     j += 1

# print(tiempo)
    #sumAGC = []
    #sumCPF = []
    #new_time = []

    #for i in range(len(time)):
    #    if t_int < time[i] < t_int + t_step:
    #        sumAGC.append(genAGC[i])
    #        sumCPF.append(genCPF[i])
    #        new_time.append(time[i])

    #P_recp = sumAGC[-1] - genAGC[0] - 292 + sumCPF[-1] - genCPF[0]
    #print('t_int = ' + str(t_int))
    #print('sumAGC = ' + str(sumAGC[-1]- 292) + ' - ' + str(genAGC[0]) + ' = ' + str(sumAGC[-1]- 292 - genAGC[0]))
    #print('sumCPF = ' + str(sumCPF[-1]) + ' - ' + str(genCPF[0]) + ' = ' + str(sumCPF[-1] - genCPF[0]))
    #print('P_recp = ' + str(P_recp))
    #P_recp_list.append(P_recp)


#sum_rocof = 0
#for t_int in range(tstop_cpf, t_final, t_step):
#     rocof=[]
#     new_time=[]

#     freq_original=[]
#     time_original=[]

#     for i in range(len(time)):
#         if t_int < time[i] < t_int + t_step:
#             freq_original.append(freq[i])
#             time_original.append(time[i])

#     percentaje = 1
#     n = int(len(time_original)*percentaje)-1  # Set the step size for calculating the slope
#     if len(time_original) >= n:
#         for j in range(n, len(time_original), n):
#             if time_original[j] - time_original[j-n] == 0:
#                 continue
#             slope = (freq_original[j] - freq_original[j-n]) / (time_original[j] - time_original[j-n])
#             rocof.append(slope)
#             new_time.append(time_original[j])

#     plt.figure()
#     plt.plot(new_time, rocof)
#     plt.title('ROCOF for t = ' + str(t_int) + ' to ' + str(t_int + t_step))
#     plt.xlabel('Time (s)')
#     plt.ylabel('ROCOF (Hz/s)')
#     plt.savefig('ROCOF_' + str(t_int) + '.png')
#     plt.close()

#     H = 1
#     sum_rocof += sum(rocof)/(2*H)

#print('Sum ROCOF = ' + str(sum_rocof))


#dif_gen_column = df['dif gen'].tolist()
#slope = df['slope'].tolist()




#print('corr1 = ' + str(df['Freq'][:2000].corr(df['dif gen'][:2000])))
#plt.figure()
#plt.scatter(df['Freq'][:2000].pct_change(),df['dif gen'][:2000].pct_change())
#plt.title('Correlation between Freq and dif gen')
#plt.savefig('correlation1.png')

#plt.figure()
#plt.scatter(df['Freq'][:2000].pct_change(),df['slope'][:2000].pct_change())
#print('corr2 = ' + str(df['Freq'][:2000].corr(df['slope'][:2000])))
#plt.title('Correlation between Freq and slope')
#plt.savefig('correlation2.png')

# Function to compare two text files
def compare_files(file1, file2):
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        f1_lines = f1.readlines()
        f2_lines = f2.readlines()
        
    differences = []
    for line_num, (line1, line2) in enumerate(zip(f1_lines, f2_lines), start=1):
        if line1 != line2:
            differences.append(f"Line {line_num}:\nFile1: {line1}\nFile2: {line2}\n")
    
    if not differences:
        print("The files are identical.")
    else:
        print("The files have the following differences:")
        for diff in differences:
            print(diff)

# Example usage
file1_path = r'C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Codigo Mauricio\Mauro\OPF_Modelo1.lp'
file2_path = r'C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\OPF_DownPower.lp'
compare_files(file1_path, file2_path)


print('Finish!')