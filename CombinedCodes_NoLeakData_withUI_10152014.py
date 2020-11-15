# -*- coding: utf-8 -*-
"""
Created on Tue Oct 07 15:54:49 2014

@author: G2BE
version 1.1.1

Change log:
10072014
- Changed indication reference from start of indication to the point with highest concentration. 
10082014
- Changed the index calculation for Wind Speed (Weibull distribution) and Intensity (50*ln(Intensity)) 
- Now both indexes are continuous.
10152014
- Changed Sigma v formula to ignore NaN data. No more missing confidence index
- Removed neighbors from CI-VI chart
10202014
- added function to export ghost indications to a csv in "Ghost Indications'.
- raised intensity minimum threshold from 15 ppb to 30 ppb
- removed leak data addition 
11032014
- changed leak file to San Jose for Picarro validation study
11042014
- changed Peak Width filter to now filter small peak widths (<1m) before multiplying with car wind angle factor (SineA)
- Now filter for very large peak Width (>40m) after applying SineA factor
11262014
- changed leak file to North Bay
- leak_near_bread gives an extra index which gives key error. 
- Most likely due to leaks close together.
- Putting another line to redefine leak_near_bread on
- added if statement to make all indications ghost if there is no leak in plat.
"""

from __future__ import division
import Tkinter
from tkFileDialog import askopenfilename
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import collections
#import scipy.stats as s
from scipy.spatial import KDTree, cKDTree
from scipy import inf

#import tkMessageBox

def open_file_handler():
    filePath = askopenfilename()
    
    return filePath
    
def getCoordinateProjectedAlongADirn(orilat,orilon,dirn,dist): #Project a coordinate from a coordinate given dirn and distance
    # dist is in km
    R = 6378.1 #Radius of the Earth
    
    lat_rad = np.radians(orilat) #Current lat point converted to radians
    lon_rad = np.radians(orilon) #Current long point converted to radians

    proj_lat = np.arcsin(np.sin(lat_rad)*np.cos(dist/R) + (
            np.cos(lat_rad)*np.sin(dist/R)*np.cos(np.radians(dirn))))

    proj_lon = lon_rad + np.arctan2(
            (np.sin(np.radians(dirn))*np.sin(dist/R)*np.cos(lat_rad)),
            (np.cos(dist/R)-np.sin(lat_rad )*np.sin(proj_lat)))

    proj_lat = np.degrees(proj_lat)
    proj_lon = np.degrees(proj_lon)
    
    return proj_lat, proj_lon

 
def getAngleUsingCosineRule(side_a,side_b,side_c):
    A = np.degrees(np.arccos(
        ((side_b)**2+(side_c)**2-(side_a)**2)/
        (2*side_b*side_c)))
    
    return A

def getScalarProductof2Vectors(startlat1, startlon1, endlat1, endlon1, 
                                startlat2, startlon2, endlat2, endlon2):
    vect1lat = endlat1-startlat1
    vect1lon = endlon1-startlon1
    vect2lat = endlat2-startlat2
    vect2lon = endlon2-startlon2
    
    dot = vect1lat * vect2lat + vect1lon * vect2lon
    
    return dot
    
def getDistanceFromLatLonInM(lat1,lon1,lat2,lon2):
    R = 6371 #Radius of the earth in km
    
    dLat = np.radians(lat2-lat1) #deg2rad below
    dLon = np.radians(lon2-lon1) 
    
    a = np.sin(dLat/2) * np.sin(dLat/2) +(
    np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * 
    np.sin(dLon/2) * np.sin(dLon/2))
    
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    d = R * c * 1000#Distance in m
    
    return d
    
def getPerpendicularCoordinates(stlat, stlon, enlat, enlon, dis):
    dx = stlat-enlat
    dy = stlon-enlon
    distreal = (dx**2+dy**2)**0.5 #not accurate formula but gives 10 meters away from stlat,stlon.
    #dist = getDistanceFromLatLonInM(stlat,stlon,enlat,enlon) # in meter
    #distreal = dist / 111111 # convert from meter to degree long/lat
    dx /= distreal
    dy /= distreal
    rightlat = stlat + dis * dy
    rightlon = stlon - dis * dx
    leftlat = stlat - dis * dy
    leftlon = stlon + dis * dx
    
    return rightlat, rightlon, leftlat, leftlon
    
def checkIfPointInParallelogram(P_lat,P_lon,A_lat,A_lon,Q_lat,Q_lon,R_lat,R_lon):#http://goo.gl/rj5I7m 
    PA_lat=A_lat-P_lat
    PA_lon=A_lon-P_lon
    PQ_lat=Q_lat-P_lat
    PQ_lon=Q_lon-P_lon
    PR_lat=R_lat-P_lat
    PR_lon=R_lon-P_lon
    detPQ_PR=PQ_lat*PR_lon-PR_lat*PQ_lon
    detPA_PQ=PA_lat*PQ_lon-PQ_lat*PA_lon
    detPA_PR=PA_lat*PR_lon-PR_lat*PA_lon
    n= -(detPA_PQ/detPQ_PR)
    m= (detPA_PR/detPQ_PR)
    #if n and m are BOTH between 0 and 1, point A is in parallelogram PQR 
    
    return n, m
    
def weib(x,n,a):   # Weibull distribution function
    return (a / n) * (x / n)**(a - 1) * np.exp(-(x / n)**a)
  

    # MAIN BODY OF CODE STARTS HERE. The code starts with the raw Picarro text file.
  
def process_file():
    filePath = open_file_handler()
    data = pd.read_csv(filePath, delim_whitespace= True)

    # Converts Epoch time to date time
    data['Date_Time'] = pd.to_datetime(data['EPOCH_TIME'],unit='s')

    # Computes the sliding average concentration by taking the last 9 and next 9 readings
    data['PeakCH4'] = pd.rolling_mean(data['CH4'], window=18).shift(-8) 
    data['PeakCH4'][:11] = pd.expanding_mean(data['CH4'])[:11] # expanding mean for the first 12 data points.
    data['PeakCH4'][-8:] = pd.rolling_mean(data['CH4'], window=9)[-8:] # Rolling mean for last few points.

    #Curve_Fit or Intensity is now the difference between point concentration and sliding average
    data['Curve_Fit']=data['CH4']-data['PeakCH4']

    #The minimum value of Curve_Fit is zero
    data['Curve_Fit'][data['Curve_Fit']<0]=0 

    #Calculate number of missing rows in car speed data
    missing= len(data.CAR_SPEED.index) - data.CAR_SPEED.count()
    pctmissing = missing / len(data.CAR_SPEED.index)
    
    #Computing Distance: cumulative distance travelled by car 
    data['Predistance'] = data['CAR_SPEED'] + data['CAR_SPEED'].shift()
    data['Predistance'] [0] = 0

    data['Predistance'][(pd.isnull(data['CAR_SPEED']))] =data['CAR_SPEED'].shift()*2
    data['Predistance'][(pd.isnull(data['CAR_SPEED'].shift()))] =data['CAR_SPEED']*2
    data['Predistance'][(pd.isnull(data['CAR_SPEED'])) & (pd.isnull(data['CAR_SPEED'].shift()))] =0

    #Distance is in (meters) = Difference in time * Speed
    data['Distance'] = ((data['EPOCH_TIME']-data['EPOCH_TIME'].shift()) * data['Predistance'] )/2
    data['Distance'][0] = 0
    data['Distance'] = data['Distance'].cumsum()

    data=data.drop(['Predistance'],1)
    
    #Compute Wind Speed
    data['Wind_Speed']=(data['WIND_N']**2+data['WIND_E']**2)**0.5
    
    #Compute Wind_Direction
    data['Wind_Dirn']=180 + np.degrees(np.arctan (data['WIND_E'] / data['WIND_N']))
    data.loc[(data['WIND_N'] >0) & (data['WIND_E'] >0)  ,'Wind_Dirn'] = data['Wind_Dirn'] - 180
    data.loc[(data['WIND_N'] >0) & (data['WIND_E'] <0)  ,'Wind_Dirn'] = data['Wind_Dirn'] + 180
    
    #Compute Car Direction
    data['Car_N']=data['GPS_ABS_LAT'] - data['GPS_ABS_LAT'].shift()
    data['Car_E']=(data['GPS_ABS_LONG']-data['GPS_ABS_LONG'].shift()) * np.cos(np.radians(data['GPS_ABS_LAT']))

    data['Car_Dirn']=180 + np.degrees(np.arctan (data['Car_E'] / data['Car_N']))
    #data.loc[(data['CAR_SPEED']<0.1) | (data['Car_N']==0), 'Car_Dirn'] = 0
    data.loc[(data['Car_N'] >0) & (data['Car_E'] >0)  ,'Car_Dirn'] = data['Car_Dirn'] - 180
    data.loc[(data['Car_N'] >0) & (data['Car_E'] <0)  ,'Car_Dirn'] = data['Car_Dirn'] + 180
    
    #Compute Car Wind Angle
    data['Car_Wind_Angle']=1
    data['min1']=abs(data['Wind_Dirn']-data['Car_Dirn'])
    data['min2']=abs(data['Wind_Dirn']-data['Car_Dirn']+360)
    data['min3']=abs(data['Wind_Dirn']-data['Car_Dirn']-360)
    data['Car_Wind_Angle'] = data[['min1','min2','min3']].min(axis=1)
    data.loc[(data['CAR_SPEED']<0.1), 'Car_Dirn'] = np.nan #Car is stopped
    data.loc[(data['CAR_SPEED']<0.1), 'Car_Wind_Angle'] = np.nan #Car is stopped
    
    #Clean up intermediate columns
    data = data.drop (['min1','min2','min3'],1)
    
    #Get mean of car wind angle column
    #data['Car_Wind_Angle'].mean()
    CWA_Skew= data['Car_Wind_Angle'].skew()
    
    #print "The original skew is:", CWA_Skew
    
    
    #Create a new column for delta time
    data['Delta_Time']=data['EPOCH_TIME']-data['EPOCH_TIME'].shift()
    
    #Compute car speed in meter per second
    data['Car_N(m/s)']=data['Car_N']*1852*60/data['Delta_Time']
    data['Car_E(m/s)']=data['Car_E']*1852*60/data['Delta_Time']
    data['Car_Speed(m/s)']=((data['Car_N(m/s)']**2)+(data['Car_E(m/s)']**2))**0.5
    
    #Compute car speed / wind speed
    Car_Speed_Mean=data['CAR_SPEED'].mean()
    Wind_Speed_Mean=data['Wind_Speed'].mean()
    Car_Wind_Ratio=Car_Speed_Mean/Wind_Speed_Mean
    
    #print 'car speed, wind speed, car/wind', Car_Speed_Mean, Wind_Speed_Mean, Car_Wind_Ratio
    
    if (abs(CWA_Skew)>0.15):
        CFRatio=CWA_Skew/2.0081
        #print"ANEMOMETER needs calibration, the Correction factor*Car-Wind-Ratio is:", CFRatio
    else:
        CFRatio=0
        #print "ANEMOMETER is in good shape"
    
    CF = CFRatio/Car_Wind_Ratio
    data['Delta_Car_N']=CF*data['Car_N(m/s)']
    data['Delta_Car_E']=CF*data['Car_E(m/s)']
    
    #CORRECTED Wind Direction
    data['Wind_N_Corr']=data['WIND_N']-data['Delta_Car_N']
    data['Wind_E_Corr']=data['WIND_E']-data['Delta_Car_E']
    
    #CORRECTED Wind Speed
    data['Wind_speed_Corr']=((data['Wind_N_Corr']**2)+(data['Wind_E_Corr']**2))**0.5
    MeanWSpeed= data['Wind_speed_Corr'].mean()    
    
    # CORRECTED Wind Direction
    data['Wind_Dirn_Corr']=180 + np.degrees(np.arctan (data['Wind_E_Corr'] / data['Wind_N_Corr']))
    data.loc[(data['Wind_N_Corr'] >0) & (data['Wind_E_Corr'] >0)  ,'Wind_Dirn_Corr'] = data['Wind_Dirn_Corr'] - 180
    data.loc[(data['Wind_N_Corr'] >0) & (data['Wind_E_Corr'] <0)  ,'Wind_Dirn_Corr'] = data['Wind_Dirn_Corr'] + 180
    
    # CORRECTED Car_Wind_Angle with NaN when CAR_SPEED is less than 0.1 km/sec
    data['Car_Wind_Angle_Corr']=2
    data['cwa1']=abs(data['Wind_Dirn_Corr']-data['Car_Dirn'])
    data['cwa2']=abs(data['Wind_Dirn_Corr']-data['Car_Dirn']+360)
    data['cwa3']=abs(data['Wind_Dirn_Corr']-data['Car_Dirn']-360)
    data['Car_Wind_Angle_Corr'] = data[['cwa1','cwa2','cwa3']].min(axis=1)
    data.loc[(data['CAR_SPEED']<0.1), 'Car_Wind_Angle_Corr'] = np.nan        
        
    #New skew in data after applying correction factor
    #CWA_Skew_Corr= data['Car_Wind_Angle_Corr'].skew()
    #print "The corrected skew is:", CWA_Skew_Corr
    

    #Filter parameters

    Min_Intensity = 0.01         # threshold for intensity in ppm
    Min_Peak_Height = 0.030 # in ppm
    #Min_Sum = 0.03          # in ppm   
    #No_Points_Checked = 6   # no. of points averaged
    Min_Peak_Sharpness = 0.5 # (Max concn-ave concn)/ave concn
    Max_Peak_Width= 40      # Max width of peaks in meters
    Min_Peak_Width= 1       # Min width of peaks in meter 
    #Dist_From_Adjacent_Peak=       #Separation from next peak, not used currently      
        
    #Create new column after applying Min. concentration filter
    data['Intensity']=0   # Intensity lists the concentration on the street
    data.loc[(data['Curve_Fit']>Min_Intensity), 'Intensity']=data['Curve_Fit']
    
    # New column: Peak detection takes the max of six readings and compare it versus the threshold
    data['Peak_Detection']=0
    data.loc[(data['Intensity']>0) & (data['Intensity'].shift()<=0) & (pd.rolling_max(data['Intensity'],window=6).shift(-5)> Min_Peak_Height) ,'Peak_Detection'] = 1
    
    data['Peak_Value']=pd.rolling_max(data['Intensity'],window=6).shift(-5)
    data['Mean_Intensity']=pd.rolling_mean(data['Intensity'],window=6).shift(-5)
    
    # Peaks are indicated by a "1" in the Peak_Detection series/array
    #Number_of_Peaks= data['Peak_Detection'][data['Peak_Detection']==1].count()
    
    # Car wind angle factor for sharpness and width filter: Sine A
    data['Sine_A'] = np.sin(np.radians(data['Car_Wind_Angle']))
    
    #Peak counter to count and mark the indications we have
    data.loc[data['Peak_Detection']>0,'Peak_Counter'] = data['Peak_Detection'].cumsum()
    
    #Sharpness filter, including the sqrt(car-wind-angle factor, sine A) * to the mean of the intensity
    #Sharpness = [Max- Mean*Sqrt(Sine A)]/ Max 
    data['Sharpness_Value']=((data['Peak_Value']-data['Mean_Intensity']*((data['Sine_A'])**0.5))/data['Peak_Value'])
    data['Sharpness_Filter'] = 0
    data.loc[(data['Peak_Detection']>0) & (data['Sharpness_Value']>Min_Peak_Sharpness), 'Sharpness_Filter'] = 1
    
    #Distance from previous peak, Caveat: the first peak has missing distance
    nonzeroid=np.nonzero(data.Peak_Counter) # return positions of non-zero elements in Peak_Counter
    Peakdistance=data.Distance[nonzeroid[0]]
    data['Peak_Dist']=Peakdistance-Peakdistance.shift()
    data.Peak_Dist.dropna()
    
    # The next section is to determine Peak width
    
    #nonzerointensity=np.nonzero(data.Intensity)
    #Peakwidth=data.Distance[nonzerointensity[0]]
    #data['Peakwidth']=Peakwidth
    data['Peak_start_locn']=np.nan#determining location of peak start
    data.loc[(data['Peak_Counter'] > 0),'Peak_start_locn']=1
    data['Peak_end_locn'] =np.nan #determining location of peak end
    data.loc[(data['Intensity'] !=0) & (data['Intensity'].shift(-1) ==0),'Peak_end_locn']=5
    data['Combined_locn']= data['Peak_start_locn'].add(data['Peak_end_locn'], fill_value=0)
    
    com = data.Combined_locn.dropna() #compressed array of start and end location
    
    data1 = pd.DataFrame(data = com) #Creating a new dataframe with markers for peak start and end positions
    data1.loc[(data1.Combined_locn == 5) & (data1.Combined_locn.shift() != 1)] = np.nan #Finding incorrect peak end locations
    data1=data1.dropna()  # Dropping rows that have incorrect peak end locations
    
    data= data.drop(['Combined_locn'],1) #Dropping column with the same name before
    data=data.join(data1)#Joining dataframes based on the index
    data['NextDistance']=data['Distance'].shift(-1) #Creating a shifted distance column for peak width calcn
    
    col_data2 = ['Combined_locn', 'NextDistance', 'Distance','Sine_A']
    data2= data.loc[pd.notnull(data['Combined_locn'])][col_data2]  #Creating new dataframe with removed rows if there is no peak start or end
    data2['Peak_Width']= np.nan   #create new column for peak width
    data2.loc[(data2['Combined_locn']==1),'Peak_Width']=data2.NextDistance.shift(-1)-data2.Distance # Calculate peak width when there is a peak start "1"
    if data2.apply(lambda x: 6 in x.values, axis=1).any():
        data2.loc[(data2['Combined_locn']==6),'Peak_Width']=data2.NextDistance-data2.Distance # Calculate peak width when there is a peak start "6"
    
    
    #Filtering with Peak Width, losing data that are bigger than Max_Peak_Width and Smaller than Min_Peak_Width
    data2.loc[data2['Peak_Width']<Min_Peak_Width, 'Peak_Width'] = np.nan #Filter out very small peak widths
    data2['Peak_Width']=((data2['Sine_A'])**0.5)*data2['Peak_Width']   #Multiply with car-wind angle factor. Square root so factor is smaller.
    data2.loc[data2['Peak_Width']>Max_Peak_Width, 'Peak_Width'] = np.nan #Filter out very large peak widths
    data2= data2.loc[pd.notnull(data2['Peak_Width'])] #removing rows without peak width information
     
    data=data.join(data2['Peak_Width'])    #Joining peak width data into main dataframe, now peak width can be found at the start of the peaks
    
    data.loc[data['Peak_Width']>0,'Peak_Detection_New'] = 1
    data.loc[data['Peak_Width']>0,'Peak_Counter_New'] = data['Peak_Detection_New'].cumsum()
    
    #NEED to add or remove IF PEAK counter
    
    #Column for average Wind Speed, rolling mean of last 30 wind speed
    WindowSpeed=30
    data.loc[(data['Peak_Counter'] != 0),'Ave_WSpeed'] = pd.rolling_mean(data['Wind_speed_Corr'].dropna(), window =WindowSpeed).reindex_like(data['Wind_speed_Corr'])
    
    #Assigning index to wind speed, 1<S<5 is 3, 0.5<S<1 is 2, S<0.5 or S>5 is 1.
    data['WSpeed_Index']= weib(data['Ave_WSpeed'],2,2)/0.3894*3
    
    #Calculation for Ave Wind Direction, using vector sum method of past 30 rows
    #Vector sum method: Convert to radians. Add sum of cosine and sum of sin. Then apply Arctan2 (ave cosine, ave sine)
    WindowDirn =30
    data['Cos_WDirn']= np.cos(np.radians(data['Wind_Dirn_Corr']))
    data['Sin_WDirn']= np.sin(np.radians(data['Wind_Dirn_Corr']))
    #np arctan2 is flipped compared to excel's (y,x) for np and (x,y) for excel 
    #Added feature to exclude NaN data while taking average direction and reduce missing directions
    data.loc[(data['Peak_Counter'] != 0),'Ave_WDirn'] = np.degrees(np.arctan2((pd.rolling_mean(data['Sin_WDirn'].dropna(), window =WindowDirn).reindex_like(data['Sin_WDirn'])),(pd.rolling_mean(data['Cos_WDirn'].dropna(), window =WindowDirn).reindex_like(data['Cos_WDirn']))))
    data.loc[(data['Ave_WDirn']<0),'Ave_WDirn'] = 360 + data['Ave_WDirn']
    
    
    #Yamartino method to calculate STDev of wind direction, Sigma_v
    data['Epsilon']= (1-((pd.rolling_mean(data['Sin_WDirn'].dropna(), window =WindowDirn).reindex_like(data['Sin_WDirn']))**2 +(
                 pd.rolling_mean(data['Cos_WDirn'].dropna(), window =WindowDirn).reindex_like(data['Cos_WDirn']))**2))**0.5
    data['Sigma_v']= np.degrees(np.arcsin(data['Epsilon'])*(1+(2/((3)**0.5)-1)*data['Epsilon']**3))
    
    
    #Assigning index to wind direction, sigma<20 is 3, 20<sigma<40 is 2, sigma >40 is 1
    data.loc[(data['Sigma_v'] <= 30), 'WDirn_Index'] = 3
    data.loc[(data['Sigma_v'] <= 60) & (data['Sigma_v'] > 20) , 'WDirn_Index'] = 2
    data.loc[(data['Sigma_v'] > 60), 'WDirn_Index'] = 1
    
    #Assigning index to maximum Intensity of an indication, 6 next cells
    data['Intens_Index']=50*np.log(data['Peak_Value']+1)
    
    #Assigning index to Plume diameter, Bigger than 45m is 0
    data.loc[(data['Peak_Width']<15),'PWidth_Index'] = 1
    data.loc[(data['Peak_Width']>15) & (data['Peak_Width']<30),'PWidth_Index'] = 0.8
    data.loc[(data['Peak_Width']>30) & (data['Peak_Width']<45),'PWidth_Index'] = 0.5
    data.loc[(data['Peak_Width']>45),'PWidth_Index'] = 0.2
    
    #Section below is to identify coordinates of max peak position
    data['Peak_Counter2']=data['Peak_Counter_New'].replace(to_replace=np.nan, method='ffill', limit = 5) #fill up next 5 rows after peak start
    idmax= data.groupby(['Peak_Counter2'])['Intensity'].idxmax() #get index of peak max
    idstart= data[pd.notnull(data['Peak_Counter_New'])].index #get index of peak starts
    
    data['max_ind_lat']= np.nan
    data['max_ind_lon']= np.nan
    
    for x in range(idmax.count()): #loop for the number of peaks in data
        startindex = idstart[x]          #position of peak start, 
        maxindex = idmax[x+1]            #position of max peak, plus 1 because index in idmax start from 1, not 0
        data['max_ind_lat'][startindex]=data['GPS_ABS_LAT'][maxindex]
        data['max_ind_lon'][startindex]=data['GPS_ABS_LONG'][maxindex]

    #Creating Indication table (data 3) (with only the necessary columns)
    col_data3=['Peak_Counter_New','Peak_Value','max_ind_lat','max_ind_lon','Ave_WSpeed','Ave_WDirn','Sigma_v','Peak_Width']
    data3= data.loc[pd.notnull(data['Peak_Width'])][col_data3]
    
    #This is to change the indication from the start of the peak to the highest point of the peak.
    data3['GPS_ABS_LAT']=data3['max_ind_lat']
    data3['GPS_ABS_LONG']=data3['max_ind_lon']
    
    max_distance = 0.00018 # Assuming lats and longs are in decimal degrees, 1 degree ~ 111km,this corresponds to 20 meters
    points = zip(data3.GPS_ABS_LAT, data3.GPS_ABS_LONG)
    tree = cKDTree(points)
    
    point_neighbors_list = [] # Put the neighbors of each point here
    
    for point in points:
        distances, indices = tree.query(point, len(points), p=2, distance_upper_bound=max_distance)
        point_neighbors = []
        for index, distance in zip(indices, distances):
            if distance == inf:
                break
            point_neighbors.append(points[index])
            point_neighbors.append(index+1)
        point_neighbors_list.append(point_neighbors)
    
    #Create a dataframe with the neighboring coordinates next to the original coordinates
    data4= pd.DataFrame(point_neighbors_list,dtype=float)
    Column_number=len(data4.columns)
    MaxNeighbor=(Column_number-2)/2
    
    
    #Replace all None with NaN
    data4=data4.fillna(value=np.nan)
    
    #Create new column for No of Neighbors
    data4['No_of_Neighbors'] = 0
    
    for i in range(2,Column_number-1,2): 
        data4.loc[pd.notnull(data4[i]),'No_of_Neighbors']= (i/2) #create new column with no. of neighbors
        temp1= pd.DataFrame(data4.loc[pd.notnull(data4[i])][i]) #create temporary dataframe to store and split tuples
        temp1[i+1]= data4.loc[pd.notnull(data4[i])][i+1]
        temp2=pd.DataFrame(temp1[i].tolist(), columns=['lat%d' % (i/2),'long%d' % (i/2)])
        temp2=temp2.join(temp1[i+1].dropna().reset_index())
        temp2['index']=temp2['index']+1
        temp2=temp2.drop(i+1,1)
        data4= data4.merge(temp2, left_on=1, right_on='index', how='outer')
        data4=data4.drop('index',1)
            
            
    data4= data4.join(pd.DataFrame(data4[0].tolist(), columns=['lat0','long0']))
    # Data 4 description
    # Original index or peak counter is column [1] 
    # Neighbor index is in column [3]. [5], [7].... depending on how many neighbors
    # Original GPS location is ['lat0'] and ['long0']
    # Neighbor's GPS coordinate are ['lat1'], ['long1'], ['lat2'], ['long2']......
    
    # Choosing the best neighbor and eliminating the rest using peak value
    # Create new data frame for peak value
    PeakVframe=pd.DataFrame(data3['Peak_Counter_New'])
    PeakVframe=PeakVframe.join(data3['Peak_Value'])
    
    # Adding Peak value to the data4 dataframe
    for i in range(2,Column_number-1,2):
        data4= data4.merge(PeakVframe, left_on=(i+1), right_on='Peak_Counter_New', how='left')
        data4['PeakValue%d' % (i/2)]=data4['Peak_Value']
        data4= data4.drop('Peak_Value',1)
        data4=data4.drop('Peak_Counter_New',1)
    
    data4 =data4.sort([1], ascending = True)
    
    #Need to identify which of peak value 1 or peakvalue 2 is higher.
    #For peaks with two or more neighbors, pick neighbor with highest peak value
    
    MaxNeighbor=int(data4['No_of_Neighbors'].max())
        
    names = []
    for i in range(1,MaxNeighbor+1):
        names.append('PeakValue%d' % i)
    
    if MaxNeighbor == 0:         
        data4['maxID']=np.nan
        data4['maxIDnumber']=np.nan
        data4['Nbr_PeakValue']=np.nan
        data3 = data3.reset_index()
        data4 = data4.reset_index()
        data3=data3.join(data4['No_of_Neighbors'])
        data3=data3.join(data4['Nbr_PeakValue']) 
    else:
        data4['maxID']= data4[names].idxmax(axis=1)
        #ID of max neighbor
        data4.loc[pd.notnull(data4['No_of_Neighbors']),'maxIDnumber']= data4['maxID'].str.lstrip('PeakValue')
        #Transfer column for No. of neighbors to data 3 (Max neighbor location & peak no.)
        data3 = data3.reset_index()
        data4 = data4.reset_index()
        data3=data3.join(data4['No_of_Neighbors'])
        #Transfer column for best Neighbor's peak value to data 3 table
        data4.loc[pd.notnull(data4['No_of_Neighbors']),'Nbr_PeakValue']=data4[names].max(axis=1)
        data3=data3.join(data4['Nbr_PeakValue'])
        
    #Transfer column for best Neighbor's location and peak no.to data3 table
    data4['Nbr_lat']=np.nan
    data4['Nbr_long']=np.nan
    data4['Nbr_PeakCounter']=np.nan    
    for row in data4.index:
        if data4['maxIDnumber'][row] > 0:
            a = int(data4['maxIDnumber'][row])
            data4['Nbr_lat'][row]= data4['lat%d' % a][row]
            data4['Nbr_long'][row]= data4['long%d' % a][row]
            data4['Nbr_PeakCounter'][row]= data4[1+2*a][row]
            
    data3=data3.join(data4['Nbr_lat'])
    data3=data3.join(data4['Nbr_long'])
    data3=data3.join(data4['Nbr_PeakCounter'])
    
    #Drop unnecessary columns
    #data4 = data4.drop(['lat0','long0','lat1','long1','lat2','long2'],1)
    
    # Do scalar product on to determine if neighbors are along wind direction.
    
    
    # Need to find a way to find latwd and longwd, a position determined from wind direction
    R = 6378.1 #Radius of the Earth
    d = 0.5 #Distance of new coordinates from original location in km\
    
    #create a new column for direction where wind is blowing TO not from
    data3.loc[data3['Ave_WDirn']>180,'WDirn_To']=data3['Ave_WDirn']-180
    data3.loc[data3['Ave_WDirn']<180,'WDirn_To']=data3['Ave_WDirn']+180
    
    #Determine a new set of coordinates based on average wind direction and 0.5 km away from original location
    data3['lat1'] = np.radians(data3['GPS_ABS_LAT']) #Current lat point converted to radians
    data3['lon1'] = np.radians(data3['GPS_ABS_LONG']) #Current long point converted to radians
    
    data3['lat2'] = np.arcsin(np.sin(data3['lat1'])*np.cos(d/R) +
        np.cos(data3['lat1'])*np.sin(d/R)*np.cos(np.radians(data3['WDirn_To'])))
    
    data3['lon2'] = data3['lon1'] + np.arctan2((np.sin(np.radians(data3['WDirn_To']))*np.sin(d/R)*np.cos(data3['lat1'])),
                (np.cos(d/R)-np.sin(data3['lat1'])*np.sin(data3['lat2'])))
    
    data3['lat2'] = np.degrees(data3['lat2'])
    data3['lon2'] = np.degrees(data3['lon2'])
    
    
    # Then we can do dot products of two arrays (latx-lat0,longy-long0),(latwd-lat0,longwd-long0)
    
    data3['Dot']= (data3['Nbr_lat']-data3['GPS_ABS_LAT'])*(data3['lat2']-data3['GPS_ABS_LAT'])+(data3['Nbr_long']-data3['GPS_ABS_LONG'])*(data3['lon2']-data3['GPS_ABS_LONG'])
    
    # Dot product is positive if vectors are pointing in similar direction, zero if they are perpendicular, and negative if the two vectors are pointing in opposite directions
    # We can check if the neighbor is downwind.
    # If dot product is +ve and Nbr is lower (lower downwind)
    # If dot product is -ve and Nbr is higher (lower downwind)
    # Question: do we want to normalize the dot product and disregard the ones close to zero.
    
    data3.loc[(data3['Dot']>0) & (data3['Peak_Value']>data3['Nbr_PeakValue']),'Lowr_DownWd']= 'yes'
    data3.loc[(data3['Dot']<0) & (data3['Peak_Value']<data3['Nbr_PeakValue']),'Lowr_DownWd']= 'yes'
    
    # Next step: Check if the indications neighboring indications can beget a downwind score.
    
    #Adding important columns into data3( Indication table)
    data3=data3.join(data[['Distance','PeakCH4','Sharpness_Value']],on='index')

    #Flagging neighbors indication with a'N'
    data3.loc[(data3['Nbr_PeakValue']>data3['Peak_Value']),'Neighbor_Indication'] ='N'
                      
    #Combining indexes into data3
    data3=data3.join(data[['WSpeed_Index','WDirn_Index','Intens_Index','PWidth_Index']],on='index')           
    
    #Assigning confidence index to Consistency (presence of neighbor) along the street.
    data3['Neighbor_Index']=1
    data3.loc[(data3['No_of_Neighbors']>0), 'Neighbor_Index']= 2 
    
    #Assigning confidence index to whether Neighbor indication is lower downwind than original indication
    # Default score is raised to 1.5 from 1 to account for our inability to check downwind/upwind when there is no neighboring indication.
    data3['LwrDownWind_Index']=1.5
    data3.loc[(data3['Lowr_DownWd']=='yes'), 'LwrDownWind_Index']= 2
    
    #Computing Total Confidence Index
    data3['Confidence_Index']= np.floor((2*data3['WSpeed_Index']*data3['WDirn_Index']*data3['Intens_Index']*data3['PWidth_Index']*data3['Neighbor_Index']*data3['LwrDownWind_Index'])/15) 
    data3.loc[(data3['Confidence_Index'] > 5),'Confidence_Index']=5
                            
    data3 = data3.drop(['lat1','lon1','WDirn_To'],1)
    data3=data3.drop(['Nbr_lat','Nbr_long','Nbr_PeakCounter','lat2','lon2','Dot','Lowr_DownWd',],1)
    
    #Import leak locations (GPS coordinates) from csv file
    leakfile = r'C:\Users\g2be\Desktop\Python trial\North Bay Pilot 5\NB Leak Locations_Complete.csv' 
    leakframe = pd.read_csv(leakfile)
        
    
    check_radius = 0.000271 # Assuming lats and longs are in decimal degrees, 1 degree ~ 111km,this corresponds to 30 meters
    leak_points = zip(leakframe.Leak_LAT, leakframe.Leak_LONG) #create KDtree list for all leaks
    leak_tree = KDTree(leak_points)                         
    bread_points = zip(data.GPS_ABS_LAT, data.GPS_ABS_LONG)  #create KDtree list for all breadcrumbs
    bread_tree = KDTree(bread_points)  
    ind_points = zip(data3.GPS_ABS_LAT, data3.GPS_ABS_LONG) #create KDtree list for all indications
    ind_tree = KDTree(ind_points) 
    
    for leakpoint in leak_points:
        dist, indc = leak_tree.query(leakpoint, len(leak_points), p=2, distance_upper_bound= max_distance)
    
    #For every leak, check if there is a nearby Picarro breadcrumb trail. 
    #Create a dataframe with leak id and position of all the breadcrumbs near a leak (<30m)
    leak_bread_frame = pd.DataFrame(columns=('leak_id','nbcrumb_LAT','nbcrumb_LONG','leak_bread_dist'))
    leak_counter = 0
    inner_counter = 0
    for leakpoint in leak_points:
        dist, indc = bread_tree.query(leakpoint, len(leak_points), p=2, distance_upper_bound= check_radius)
        for index, distance in zip(indc, dist):
            if distance == inf:
                break
            leak_bread_frame.loc[inner_counter] =[leak_counter, data.GPS_ABS_LAT[index], data.GPS_ABS_LONG[index], distance]
            inner_counter += 1
        leak_counter += 1
    
    #Based on leak id add leak location to column and 
    
    #Using set to identify unique numbers in leak ids that are near breadcrumbs
    leak_near_bread = set(leak_bread_frame['leak_id']) 
    
    #Find leaks that are near indications 
    leak_ind_frame = pd.DataFrame(columns=('leak_id','nind_LAT','nind_LONG','leak_ind_dist'))
    inner_counter = 0
    leak_counter = 0
    for leakpoint in leak_points:
        dist, indc = ind_tree.query(leakpoint, len(leak_points), p=2, distance_upper_bound= check_radius)
        for index, distance in zip(indc, dist):
            if distance == inf:
                break
            leak_ind_frame.loc[inner_counter] =[leak_counter, data3.GPS_ABS_LAT[index], data3.GPS_ABS_LONG[index], distance]
            inner_counter += 1
        leak_counter += 1
    
    #For each pair of leak to indication, assign a Validation index 
    #Using set to identify unique numbers in leak ids that are near indications
    leak_near_ind = set(leak_ind_frame['leak_id']) #Check if ind are in wind sector
    #leak_wo_ind= leak_near_bread - leak_near_ind #Check scalar prd of leak-wind vs. leak-crumb 
    
    #For leaks without indication, find breadcrumb that has the closest direction to the wind.
    startLat=np.radians(list(leak_bread_frame['nbcrumb_LAT'])) #breadcrumb locn
    startLong= np.radians(list(leak_bread_frame['nbcrumb_LONG']))
    endLong= np.radians(list(leakframe['Leak_LONG'][leak_bread_frame['leak_id']]))  #leak locn
    endLat= np.radians(list(leakframe['Leak_LAT'][leak_bread_frame['leak_id']]))
                       
    dLong = endLong - startLong
    dPhi = np.log(np.tan(endLat/2.0+np.pi/4.0)/np.tan(startLat/2.0+np.pi/4.0))
    
    dLongn = []
    for x in dLong:
        if abs(x) > np.pi:
            if x > 0.0:
                x = -(2.0 * np.pi - x)
                dLongn.append(x)
            else:
                x = (2.0 * np.pi + x)
                dLongn.append(x)
            bearing = (np.degrees(np.arctan2(dLong, dPhi)) + 360.0) % 360.0;
    
    # bearing is a list with all the breadcrumb direction to the leak
    #Now need to get wind direction at every breadcrumb. Need to import this from 'data' dataframe [Wind_Dirn_Corr']
    # use the bcrumb lat posn to get wind direction [nbcrumb_LAT] in leak_bread_frame with [GPS_ABS_LAT]
    
    
    
    #Cr8 new dataframe to store wind dirn info for each breadcrumb before merging to leak_bread_frame with matching GPS locn
    bread_wind_frame = data[['Wind_Dirn_Corr','GPS_ABS_LAT']]
    
    leak_bread_frame = leak_bread_frame.merge(bread_wind_frame, how='inner',left_on='nbcrumb_LAT', right_on='GPS_ABS_LAT') #add Wind_Dirn to leak_bread_frame
    leak_bread_frame = leak_bread_frame[pd.notnull(leak_bread_frame['Wind_Dirn_Corr'])] # Remove breadcrumbs with NaN Wind Dirn
    leak_bread_frame = leak_bread_frame.sort(column='leak_id', ascending=True) #sort leak_bread_frame based on ascending leak_id
    
    #Insert position of leaks into leak_bread_frame to prepare for scalar product
    leak_posn_frame = leakframe[['Leak_LAT','Leak_LONG']] 
    leak_posn_frame['leak_id'] = leak_posn_frame.index
    leak_bread_frame = leak_bread_frame.merge(leak_posn_frame, how='inner', left_on='leak_id', right_on='leak_id')
    
    # Project a coordinate from leak to where the wind is blowing FROM, 0.5 km away in degree
    leak_bread_frame['wdproj_lat'],leak_bread_frame['wdproj_lon'] = getCoordinateProjectedAlongADirn(
                                                leakframe['Leak_LAT'],leakframe['Leak_LONG'],
                                                leak_bread_frame['Wind_Dirn_Corr'],0.5)
                                            
    leak_bread_frame['Dot_Prd']= getScalarProductof2Vectors(
                        leak_bread_frame['nbcrumb_LAT'], leak_bread_frame['nbcrumb_LONG'],
                        leak_bread_frame['wdproj_lat'], leak_bread_frame['wdproj_lon'],
                        leak_bread_frame['nbcrumb_LAT'], leak_bread_frame['nbcrumb_LONG'],
                        leak_bread_frame['Leak_LAT'], leak_bread_frame['Leak_LONG'])
                    
    #leak_bread_frame = leak_bread_frame.drop(['lk_lat_rad','lk_lon_rad'],1)
              
    #Find max dot product for each set of leak-bcrumb
    lk_maxdotprd = leak_bread_frame.groupby(['leak_id'])['Dot_Prd'].max()
    
    #Add column at leak dataframe for validation index, default value is Out of Bounds
    leakframe['Val_ID']='OOB'
    
    #Assign validation index, if max dot prd is +ve, val id =2, otherwise it's a 1
    leak_near_bread = set(leak_bread_frame['leak_id']) #recalculate leak_near_bread to eliminate extras
    for x in leak_near_bread:
        if lk_maxdotprd.loc[x] > 0:  #Changed to .loc to avoid key error , Nov 26, 2014
            leakframe['Val_ID'][x] = 2
        else:
            leakframe['Val_ID'][x] = 1    
     
     #Include an if statement here to account for cases when there is no leak within 30m of an indication
    if leak_ind_frame.empty:
        lk_val_id = []
        lk_con_id = []
    
    else:
        # Importing Wind Direction data into leak_indication_frame
        ind_wind_frame = data3 [['Ave_WDirn','Sigma_v','GPS_ABS_LAT','Confidence_Index','Peak_Counter_New','Neighbor_Indication']]
        leak_ind_frame = leak_ind_frame.merge(ind_wind_frame, how='inner',left_on='nind_LAT', right_on='GPS_ABS_LAT') #add Wind_Dirn to leak_ind_frame
        leak_ind_frame = leak_ind_frame[pd.isnull(leak_ind_frame['Neighbor_Indication'])] # Remove rows with "N" = Neighboring Indication
        leak_ind_frame = leak_ind_frame[pd.notnull(leak_ind_frame['Ave_WDirn'])] # Remove ind with NaN Wind Dirn
        leak_ind_frame = leak_ind_frame.sort(column='leak_id', ascending=True) #sort leak_ind_frame based on ascending leak_id
        #Insert position of leaks into leak_ind_frame
        leak_ind_frame = leak_ind_frame.merge(leak_posn_frame, how='inner', left_on='leak_id', right_on='leak_id')
        
        
        # Define wedge boundaries (wind sector) for each leak with indications
        # Define vertices of triangle TRIANGLE METHOD
        
        leak_ind_frame['nind_LAT']=np.float64(leak_ind_frame['nind_LAT']) # convert nind coordinates to float
        leak_ind_frame['nind_LONG']=np.float64(leak_ind_frame['nind_LONG'])
        
        leak_ind_frame['base_lat'], leak_ind_frame['base_lon'] = getCoordinateProjectedAlongADirn(
                                        leak_ind_frame['nind_LAT'],leak_ind_frame['nind_LONG'],
                                        leak_ind_frame['Ave_WDirn'],0.03)
                                        
        #Find Direction of Vertices of wind sector (triangle)
        leak_ind_frame['min_proj_dirn']= leak_ind_frame['Ave_WDirn']-leak_ind_frame['Sigma_v']
        leak_ind_frame['min_proj_dirn'] = np.where(leak_ind_frame['min_proj_dirn'] < 0, 360 + leak_ind_frame['min_proj_dirn'],leak_ind_frame['min_proj_dirn'])
        leak_ind_frame['max_proj_dirn']= leak_ind_frame['Ave_WDirn']+leak_ind_frame['Sigma_v']
        leak_ind_frame.loc[(leak_ind_frame['max_proj_dirn'] > 360), 'max_proj_dirn'] = leak_ind_frame['max_proj_dirn'] - 360
        
        # Determine vertices of wind sector in degrees
            
        leak_ind_frame['wdsect_vert2_lat'], leak_ind_frame['wdsect_vert2_lon'] = getCoordinateProjectedAlongADirn(
                                                leak_ind_frame['nind_LAT'],leak_ind_frame['nind_LONG'],
                                                leak_ind_frame['min_proj_dirn'],0.03)
                                                
        leak_ind_frame['wdsect_vert3_lat'], leak_ind_frame['wdsect_vert3_lon'] = getCoordinateProjectedAlongADirn(
                                                leak_ind_frame['nind_LAT'],leak_ind_frame['nind_LONG'],
                                                leak_ind_frame['max_proj_dirn'],0.03) 
                                            
        #leak_ind_frame['wdproj_lat'], leak_ind_frame['wdproj_lon'] = getCoordinateProjectedAlongADirn(
        # leak_ind_frame['nind_LAT'],leak_ind_frame['nind_LONG'],
        #leak_ind_frame['Ave_WDirn'],0.09)         
                                            
                                            
        # Check value of 2x sigma_v
        # If sigma_v is between 90 and 120, just take the scalar product of ind-leak to leak-wind, +ve is 5, -ve is 1.Check distance too.
        Threshold_circle_sector = 120 #if stdev wind is between this and 90, round down stdev to 90.
        
        # Dot prd of indn-leak, ind-wind   
        leak_ind_frame['Dot_IL_IW']= getScalarProductof2Vectors(
                        leak_ind_frame['nind_LAT'], leak_ind_frame['nind_LONG'], 
                        leak_ind_frame['Leak_LAT'], leak_ind_frame['Leak_LONG'], 
                        leak_ind_frame['nind_LAT'], leak_ind_frame['nind_LONG'], 
                        leak_ind_frame['base_lat'], leak_ind_frame['base_lon'])
                    
                    
        # Leaks downwind and when wind sector is small, val =1
        leak_ind_frame.loc[(0> leak_ind_frame['Dot_IL_IW']) & (leak_ind_frame['Sigma_v'] <90), 'Val_ID'] = 1
        
        # When wind sector is larger than 120 deg. threshold, val = 5
        leak_ind_frame.loc[(leak_ind_frame['Sigma_v']>Threshold_circle_sector),'Val_ID'] = 5
        
        # When wind sector is between 90 and 120, and leak is downwind, val =1 
        leak_ind_frame.loc[(leak_ind_frame['Sigma_v']>90) & (leak_ind_frame['Sigma_v']<Threshold_circle_sector) & (
                         0 > leak_ind_frame['Dot_IL_IW']),'Val_ID']= 1
        # when wind sector is between 90 and 120, and leak is upwind, val= 5
        leak_ind_frame.loc[(leak_ind_frame['Sigma_v']>90) & (leak_ind_frame['Sigma_v']<Threshold_circle_sector) & (
                         leak_ind_frame['Dot_IL_IW']>0),'Val_ID']= 5
                         
        #Cases when sigma_v is less than 90 degrees\
        # Cosine rule to calculate angle between leaks and vertices
        leak_ind_frame['Len_ind_lk'] = getDistanceFromLatLonInM(leak_ind_frame['nind_LAT'],leak_ind_frame['nind_LONG'],leak_ind_frame['Leak_LAT'],leak_ind_frame['Leak_LONG'])
        leak_ind_frame['Len_ind_v2'] = getDistanceFromLatLonInM(leak_ind_frame['nind_LAT'],leak_ind_frame['nind_LONG'],leak_ind_frame['wdsect_vert2_lat'],leak_ind_frame['wdsect_vert2_lon'])
        leak_ind_frame['Len_v2_lk'] = getDistanceFromLatLonInM(leak_ind_frame['Leak_LAT'],leak_ind_frame['Leak_LONG'],leak_ind_frame['wdsect_vert2_lat'],leak_ind_frame['wdsect_vert2_lon'])
        leak_ind_frame['Ang_v2_id_lk'] = getAngleUsingCosineRule(leak_ind_frame['Len_v2_lk'],leak_ind_frame['Len_ind_lk'],leak_ind_frame['Len_ind_v2'])
    
    
        leak_ind_frame['Len_ind_v3'] = getDistanceFromLatLonInM(leak_ind_frame['nind_LAT'],leak_ind_frame['nind_LONG'],leak_ind_frame['wdsect_vert3_lat'],leak_ind_frame['wdsect_vert3_lon'])
        leak_ind_frame['Len_v3_lk'] = getDistanceFromLatLonInM(leak_ind_frame['Leak_LAT'],leak_ind_frame['Leak_LONG'],leak_ind_frame['wdsect_vert3_lat'],leak_ind_frame['wdsect_vert3_lon'])
        leak_ind_frame['Ang_v3_id_lk'] = getAngleUsingCosineRule(leak_ind_frame['Len_v3_lk'],leak_ind_frame['Len_ind_lk'],leak_ind_frame['Len_ind_v3'])

        
        # leak is in sector if angles between leak-ind-vert are less than total wind variation
        leak_ind_frame.loc[(leak_ind_frame['Ang_v2_id_lk']< (2*leak_ind_frame['Sigma_v']))& (
                        leak_ind_frame['Ang_v3_id_lk'] < (2*leak_ind_frame['Sigma_v'])) & (
                        leak_ind_frame['Sigma_v']<90),'Val_ID'] = 5
                        
        leak_ind_frame = leak_ind_frame.drop(['Ang_v2_id_lk','Ang_v3_id_lk','Len_ind_lk', 'Len_ind_v2', 'Len_v2_lk', 'Len_ind_v3', 'Len_v3_lk','GPS_ABS_LAT'],1)
        
        # Count leaks that are 10m away from sector = Val ID 4
        
        # Define base position 10 meters away from center of circle on both sides. Base of parallelogram
        N =  0.00009033# distance from indication, 10 meters
    
        leak_ind_frame['rightlat'], leak_ind_frame['rightlon'], leak_ind_frame['leftlat'], leak_ind_frame['leftlon'] = getPerpendicularCoordinates(
        leak_ind_frame['nind_LAT'], leak_ind_frame['nind_LONG'], leak_ind_frame['base_lat'], leak_ind_frame['base_lon'], N)
         
        leak_ind_frame['n_l'], leak_ind_frame['m_l']= checkIfPointInParallelogram(leak_ind_frame['nind_LAT'], leak_ind_frame['nind_LONG'],
                                leak_ind_frame['Leak_LAT'], leak_ind_frame['Leak_LONG'], 
                                leak_ind_frame['wdsect_vert2_lat'], leak_ind_frame['wdsect_vert2_lon'],
                                leak_ind_frame['leftlat'], leak_ind_frame['leftlon'])
                                    
        leak_ind_frame['n_r'], leak_ind_frame['m_r']= checkIfPointInParallelogram(leak_ind_frame['nind_LAT'], leak_ind_frame['nind_LONG'],
                                leak_ind_frame['Leak_LAT'], leak_ind_frame['Leak_LONG'], 
                                leak_ind_frame['wdsect_vert3_lat'], leak_ind_frame['wdsect_vert3_lon'],
                                leak_ind_frame['rightlat'], leak_ind_frame['rightlon'])
                                
        leak_ind_frame.loc[(0 <= leak_ind_frame['n_r']) & (leak_ind_frame['n_r'] <= 1) & (0 <= leak_ind_frame['m_r']) & (leak_ind_frame['m_r']<=1), 'rcheck'] = 1
        leak_ind_frame.loc[(0 <= leak_ind_frame['n_l']) & (leak_ind_frame['n_l'] <= 1) & (0 <= leak_ind_frame['m_l']) & (leak_ind_frame['m_l']<=1), 'lcheck'] = 1
        
        leak_ind_frame.loc[((leak_ind_frame['rcheck'] == 1) | (leak_ind_frame['lcheck'] == 1)) & (
                            leak_ind_frame['Sigma_v']<90), 'Val_ID'] = 4
                        
        leak_ind_frame.drop(['n_r','n_l','m_r','m_l'],1)
        
        #If a leak is within 5 meters of indication, Val ID is a 5
        leak_ind_frame.loc[(leak_ind_frame['leak_ind_dist']<0.000045),'Val_ID'] = 5
        
        
        #Anything that is not a 5,4,1 is a 3
        leak_ind_frame.loc[(leak_ind_frame['Val_ID'] != 5) & (leak_ind_frame['Val_ID'] !=4) & (leak_ind_frame['Val_ID'] !=1),'Val_ID']=3
        
        
        #If we do not take the highest Val ID and the highest Con ID
        lk_val_id=leak_ind_frame['Val_ID']
        lk_con_id=leak_ind_frame['Confidence_Index']
        
        

    #THIS SECTION IS TO PLOT GRAPHS and CHARTS


    fig1 = plt.figure(1)
    gs1=gridspec.GridSpec(4,3)
    
    # Creating the GridSpec to arrange multiple plots in a figure
    
    
    #Car Direction distribution histogram, plot 5
    ax5 = fig1.add_subplot(gs1[1, -1])
    x, bins, p = plt.hist(data['Car_Dirn'], bins=6, range= (0,360), facecolor='green',alpha = 0.7, align='left',rwidth=0.7)
    plt.xticks(np.arange(0,360,60))              # Setting the ticks on x axis
    plt.ylim(0,0.4)                             # Setting the limit of y axis
    #plt.xlabel('Car Direction (degrees)')
    plt.ylabel('Fraction')
    plt.title(r'Car Direction')
    for item in p:
        item.set_height(item.get_height()/sum(x)) # Normalizing the data so the sum of the bar height is =1
    
    
    #Wind Speed distribution histogram, plot 4
    ax4 = fig1.add_subplot(gs1[1,-2])
    x, bins, p = plt.hist(data['Wind_Speed'], bins=8, range= (0,4),alpha=0.4, facecolor='blue',align='left',rwidth=0.55,label = 'Original')
    y, ybins, py = plt.hist(data['Wind_speed_Corr'], bins=8, range= (0,4),alpha=0.4, facecolor='yellow',align='left',rwidth=0.7,label = 'Corrected')
    plt.xticks(np.arange(0,4,0.5))
    plt.ylim(0,0.55)
    plt.annotate ('Ave= %.2f m/s' % MeanWSpeed, xy=(0.55,0.626),xycoords='figure fraction', arrowprops=None, fontsize =11)
    #plt.xlabel('Wind Speed (m/s)')
    plt.ylabel('Fraction')
    plt.title(r'Wind Speed (m/s)')
    for item in p:
        item.set_height(item.get_height()/sum(x))
    for item in py:
        item.set_height(item.get_height()/sum(y))
    plt.legend(loc=1,prop={'size':10})
    
    
    #Car-Wind-Angle distribution histogram, plot 2
    ax2 = fig1.add_subplot(gs1[0, 2])
    x, bins, p = plt.hist(data['Car_Wind_Angle'], bins=6, range= (0,180),alpha=0.4, facecolor='red',align='mid',rwidth=0.55, label='Original')
    y, ybins, py = plt.hist(data['Car_Wind_Angle_Corr'], bins=6, range= (0,180),alpha=0.4, facecolor='blue',align='mid',rwidth=0.7,label = 'Corrected')
    plt.xticks(np.arange(0,180,30))
    plt.ylim(0,0.4)
    plt.annotate('Correction factor = %.2f' % CF, xy=(0.27,0.9),xycoords='figure fraction', arrowprops=None,fontsize=13)
    plt.ylabel('Fraction')
    plt.title(r'Car-Wind Angle')
    for item in p:
        item.set_height(item.get_height()/sum(x))
    for item in py:
        item.set_height(item.get_height()/sum(y))
    plt.legend(loc=1,prop={'size':10})
    
        
    #CH4 Concentration profile line plot, plot 7
    ax7 = fig1.add_subplot(gs1[2, 1:])
    plt.plot(data['Distance'], data['CH4'])
    plt.ylim(1.5, 4.5) # Before was plt.ylim(plt.ylim()[0], 4.5)
    plt.ylabel('CH4 Concentration (ppm)')
    plt.show()
    
    #Intensity profile over distance, line plot, plot 9
    ax9 = fig1.add_subplot(gs1[-1, 1:])
    plt.plot(data['Distance'], data['Curve_Fit'], color = 'brown')
    plt.ylim(-0.05, 0.6)
    plt.xlabel('Distance (m)')
    plt.ylabel('Intensity (ppm)')
    plt.show()
    
    
    #Draw Table with Concentration and Curve Fit data, plot 1
    ax1 = fig1.add_subplot(gs1[0, :-1])
    col_labels=['Concentration','(ppm)','Intensity','(ppm)']
    table_vals= [                                                 # Manually insert table values
        ['Min',round(data['CH4'].min(),3),'Threshold','>0.01 ppm'],
     ['Max',round(data['CH4'].max(),3),'Min Peak',round(data['Curve_Fit'].min(), 3)],
     ['Ave',round(data['CH4'].mean(),3),'Max Peak',round(data['Curve_Fit'].max(),3)],
     ['Median',round(data['CH4'].median(),3),'Above threshold','%d counts' % round(data.Curve_Fit[data.Curve_Fit>0.01].count(),3)],
     ['Ave-Min',round(data['CH4'].mean()-data['CH4'].min(),3),'Ratio','%.2f %%' % (100*(data.Curve_Fit[data.Curve_Fit>0.01].count()/data.Curve_Fit.count()))],
     ['Missing data','%d rows' % missing, '% missing data', '%.2f %%' % (100*pctmissing)]]
    
    hcell, wcell = 0.5, 1        # Set column width and row height

    the_table = plt.table(cellText=table_vals,
                      colLabels=col_labels,
        loc = 'center', colColours=['skyblue','skyblue','orange','orange'])     # Plot table
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(11)
    the_table.scale(0.9, 0.9)
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.axis('off')     #Hide x and y axis
    path, file_name = os.path.split(filePath) #splits path and file name for title 
    plt.annotate ('Condition for '+ file_name[:-29], xy=(0.15,0.95),xycoords='figure fraction', arrowprops=None, fontsize =17, color = 'crimson')

    # Draw Wind Direction Distribution Radar Plot, Plot 3
    # First count the number of instances for every 30 degree sector
    totalwinddirncorr = data.Wind_Dirn_Corr.count()
    ccount1 = data.Wind_Dirn_Corr[data.Wind_Dirn_Corr<=30].count()/totalwinddirncorr
    ccount2 = data.Wind_Dirn_Corr[(data.Wind_Dirn_Corr>30) & (data.Wind_Dirn_Corr<=60)].count()/totalwinddirncorr
    ccount3 = data.Wind_Dirn_Corr[(data.Wind_Dirn_Corr>60) & (data.Wind_Dirn_Corr<=90)].count()/totalwinddirncorr
    ccount4 = data.Wind_Dirn_Corr[(data.Wind_Dirn_Corr>90) & (data.Wind_Dirn_Corr<=120)].count()/totalwinddirncorr
    ccount5 = data.Wind_Dirn_Corr[(data.Wind_Dirn_Corr>120) & (data.Wind_Dirn_Corr<=150)].count()/totalwinddirncorr
    ccount6 = data.Wind_Dirn_Corr[(data.Wind_Dirn_Corr>150) & (data.Wind_Dirn_Corr<=180)].count()/totalwinddirncorr
    ccount7 = data.Wind_Dirn_Corr[(data.Wind_Dirn_Corr>180) & (data.Wind_Dirn_Corr<=210)].count()/totalwinddirncorr
    ccount8 = data.Wind_Dirn_Corr[(data.Wind_Dirn_Corr>210) & (data.Wind_Dirn_Corr<=240)].count()/totalwinddirncorr
    ccount9 = data.Wind_Dirn_Corr[(data.Wind_Dirn_Corr>240) & (data.Wind_Dirn_Corr<=270)].count()/totalwinddirncorr
    ccount10 = data.Wind_Dirn_Corr[(data.Wind_Dirn_Corr>270) & (data.Wind_Dirn_Corr<=300)].count()/totalwinddirncorr
    ccount11 = data.Wind_Dirn_Corr[(data.Wind_Dirn_Corr>300) & (data.Wind_Dirn_Corr<=330)].count()/totalwinddirncorr
    ccount12 = data.Wind_Dirn_Corr[(data.Wind_Dirn_Corr>330) & (data.Wind_Dirn_Corr<=360)].count()/totalwinddirncorr
    
    windarray = [ccount1, ccount2, ccount3, ccount4, ccount5, ccount6, ccount7, ccount8, ccount9, ccount10, ccount11, ccount12, ccount1]
    theta = np.arange(0, 360, 30, dtype=float) * np.pi / 180.0
    thetanew = np.append(theta,[0])
    
    ax3 = fig1.add_subplot(gs1[1,0], polar = True) # plotting the radar plot
    
    plt.plot(thetanew, windarray, color= 'red')
    plt.fill(thetanew, windarray, facecolor= 'red', alpha=0.25)
    ax3.set_theta_zero_location('N')             #set 0 at north
    ax3.set_theta_direction(-1)                  # move in clockwise direction
    ax3.set_ylim(0,0.3)
    ax3.set_yticks(np.arange(0,0.3,0.1))
    plt.title(r'Wind Direction (corrected)', size=11)
        
    #IF statement here to account for survey runs with no within 30 m of an indication
    if leak_ind_frame.empty: 
        ax8 = fig1.add_subplot(gs1[2:, 0])
        plt.axis([0, 1, 0, 1])
        t = "There is no leak near"
        t2= "an indication in this run!"
        plt.text(0.5, 0.55, t, fontsize=18, family='serif', style='italic', ha='center', va='center')
        plt.text(0.5, 0.45, t2, fontsize=18, family='serif', style='italic', ha='center', va='center')
    
    else:     
        #Plot histogram of confidence indexes, plot 6
        
        ax6 = fig1.add_subplot(gs1[2, 0])
        x, bins, p = plt.hist(lk_con_id, bins=6, range= (-0.5,5.5),alpha=0.4, facecolor='red',align='mid',rwidth=0.55)
        plt.xticks(np.arange(0,6,1.0))
        plt.xlim(-0.5, 5.5)
        plt.ylim(0.05,0.55)
        plt.ylabel('% of Indications')
        plt.xlabel('Confidence Index', size = 9)
        plt.title(r'Confidence Index distribution', size =11)
        for item in p:
            item.set_height(item.get_height()/sum(x))
            
            
        #Plot bubble plot of Val ID vs. Con ID, plot 8
        scatterdata = zip(lk_con_id,lk_val_id)
        count = collections.Counter(scatterdata) #count no. of total sets (no repeats)
           
        points = count.keys() #remove repeats in data
        conx, valy = zip(*points) #split up the x and y after removing repeats
        size = np.array(count.values()) #gives count of each bubble
        sizes = 200* np.array(count.values()) # determine size of bubbles
        
        
        ax8= fig1.add_subplot(gs1[3, 0])
        plt.scatter(conx,valy, marker ='o', c='royalblue', s=sizes, alpha =0.85)
        for label, x, y in zip (size, conx, valy):
            plt.annotate(
                label, 
                xy = (x, y),color = 'white',
                textcoords='data', ha='left', size=9)
            plt.xlabel ('Confidence Index')
            plt.ylabel ('Validation Index')
            plt.ylim(-0.5, 5.5)
        plt.xlim(-0.5, 5.5)
        
        ax8.set_xticks(np.arange(-0.5, 5.5, 1),minor = True)
        ax8.set_yticks(np.arange(-0.5, 5.5, 1),minor = True)
        ax8.grid(which='minor')
    
    
    # Change font sizes of all axes and titles
    font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 11}
    plt.rc('font', **font)
    #gs1.tight_layout(fig1)    
    
    plt.show()

    path_d = 'C:\Users\g2be\Desktop\Python trial\Output' # Assign path directory for csv output.
    data3.to_csv(os.path.join(path_d, file_name + '.csv'))   
    
    # This part stores all ghost indications into a csv file and includes the confidence index information
    if leak_ind_frame.empty:
        ghost_frame = data3

    else:
        ghostidx = data3['Peak_Counter_New'].isin(leak_ind_frame['Peak_Counter_New'])
        ghost_frame = data3[~ghostidx]

    path_g = 'C:\Users\g2be\Desktop\Python trial\Ghost Indication' # Assign path directory for csv output.
    ghost_frame.to_csv(os.path.join(path_g, file_name + 'GHOST.csv'))
    
    
    return data3
    return file_name 

#This part is to generate the pop-up window user interface to select files. 
main_window = Tkinter.Tk()
        
open_file = Tkinter.Button(main_window, command=process_file, padx=150,pady=10, text="SELECT .txt FILE to process").grid(row=0)
Tkinter.Label(main_window, text="Please change output path directory in code manually").grid(row=1)

main_window.mainloop() 
