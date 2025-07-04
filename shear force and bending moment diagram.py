import math
import numpy as np
import matplotlib.pyplot as plt 

pointLoads=np.array([[]]) #point forces[location,xMagnitude,yMagnitude]
pointMoments=np.array([[]])#moments[location, magnitude]
UDL=np.array([[]])#UDL[start pos,end pos,mag]
UVL=np.array([[]])#UVL[start pos,end pos, star mag,end mag]

#input span data and force data
span=6
A=0 #left support distance
B=6 #right support distance
#pointLoads=np.array([[4,0,-60]])
#pointMoments=np.array([[7,1200]])
UDL=np.array([[0,2,-50],[4,6,-20]])
#UVL=np.array([[0,6,0,-270]])


#initialization and defualt values
div=10000
delta=span/div
X=np.arange(0,span+delta,delta)
nPL=len(pointLoads[0])
nPM=len(pointMoments[0])
nUDL=len(UDL[0])
nUVL=len(UVL[0])
reactions=np.array([0.0,0,0])
shearForce=np.empty([0,len(X)])
bendingMoment=np.empty([0,len(X)])


#define a function to calculate reactions due to point loads
def reactions_PL(n):
    xp=pointLoads[n,0]
    fx=pointLoads[n,1]
    fy=pointLoads[n,2]
    
    la_p=A-xp
    mp=fy*la_p
    la_vb=B-A
    vb=mp/la_vb
    va=-fy-vb
    ha=-fx
    return va,vb,ha
#Define a function to calculate the reaction due to the point moment
def reactions_PM(n):
    xm=pointMoments[n,0]
    m=pointMoments[n,1]
    la_vb=B-A
    Vb=m/la_vb
    Va=-Vb
    return Va,Vb

#define a function to calculate reactions due to UDL
def reactions_UDL(n):
    xstart=UDL[n,0]
    xend=UDL[n,1]
    fy=UDL[n,2]
    
    fy_resultant=fy*(xend-xstart)
    x_resultant=xstart+0.5*(xend-xstart)
    la_p=A-x_resultant
    mp=fy_resultant*la_p
    la_vb=B-A
    vb=mp/la_vb
    va=-fy_resultant-vb
    return va,vb

#define a function to calculate reactions due to UVL
def reactions_UVL(n):
    xstart=UVL[n,0]
    xend=UVL[n,1]
    fy_start=UVL[n,2]
    fy_end=UVL[n,3]
    if abs(fy_start>0):
        fy_resultant=0.5*fy_start*(xend-xstart)
        x_resultant=xstart+(1/3)*(xend-xstart)
    else:
        fy_resultant=0.5*fy_end*(xend-xstart)
        x_resultant=xstart+(2/3)*(xend-xstart)
    la_p=A-x_resultant
    mp=fy_resultant*la_p
    la_vb=B-A
    vb=mp/la_vb
    va=-fy_resultant-vb
    return va,vb
        
   
    
     

#calculating reaction due to point loads
PL_record=np.empty([0,3])
if (nPL>0): 
    for n, p in enumerate (pointLoads):
        va,vb,ha=reactions_PL(n)
        PL_record=np.append(PL_record,[np.array([va,ha,vb])],axis=0)
        #add reactions to record
        reactions[0]=reactions[0]+va
        reactions[1]=reactions[1]+ha
        reactions[2]=reactions[2]+vb
                
#calculating reaction due to point Moments
PM_record=np.empty([0,2])
if(nPM>0):
    for n,p in enumerate(pointMoments):
        va,vb=reactions_PM(n)
        PM_record=np.append(PM_record,[np.array([va,vb])],axis=0)
         #add reactions to record
        reactions[0]=reactions[0]+va
        reactions[2]=reactions[2]+vb
        
       
#calculating reaction due to UDL
UDL_record=np.empty([0,2])
if(nUDL>0):
    for n,p in enumerate(UDL):
        va,vb=reactions_UDL(n)
        UDL_record=np.append(UDL_record,[np.array([va,vb])],axis=0)
         #add reactions to record
        reactions[0]=reactions[0]+va
        reactions[2]=reactions[2]+vb

#calculating reaction due to UVL
UVL_record=np.empty([0,2])
if(nUVL>0):
    for n,p in enumerate(UVL):
        va,vb=reactions_UVL(n)
        UVL_record=np.append(UVL_record,[np.array([va,vb])],axis=0)
         #add reactions to record
        reactions[0]=reactions[0]+va
        reactions[2]=reactions[2]+vb
        
        






        
        
print("The vertical reaction at A is ",reactions[0],"kN")   
print("The vertical reaction at B is {one} kN".format(one=reactions[2]))
#print("The horizontal reaction at A is {one} kN".format(one=reactions[1]))   

#Shear and moment calculation due to point loads
def shear_moment_PL(n):
    xp=pointLoads[n,0]
    fy=pointLoads[n,2]
    Va=PL_record[n,0]
    Vb=PL_record[n,2]
    Shear=np.zeros(len(X))
    Moment=np.zeros(len(X))
    
    for i,x in enumerate(X):
        
        shear=0
        moment=0
        if x>A:
            shear=shear+Va
            moment=moment-Va*(x-A)
            
        if x>B:
            shear=shear+Vb
            moment=moment-Vb*(x-B)
        
        if x>xp:
            shear=shear+fy
            moment=moment-fy*(x-xp)
            
        Shear[i]=shear
        Moment[i]=moment
        
    return Shear , Moment  


#Shear and moment calculation due to point Moments
def shear_moment_PM(n):
    xm=pointMoments[n,0]
    m=pointMoments[n,1]
    Va=PM_record[n,0]
    Vb=PM_record[n,1]
    Shear=np.zeros(len(X))
    Moment=np.zeros(len(X))
    
    for i,x in enumerate(X):
        
        shear=0
        moment=0
        if x>A:
            shear=shear+Va
            moment=moment-Va*(x-A)
            
        if x>B:
            shear=shear+Vb
            moment=moment-Vb*(x-B)
        
        if x>xm:
            
            moment=moment-m
            
        Shear[i]=shear
        Moment[i]=moment
        
    return Shear , Moment  

#Shear and moment calculation due to UDL
def shear_moment_UDL(n):
    xstart=UDL[n,0]
    xend=UDL[n,1]
    fy=UDL[n,2]
    Va=UDL_record[n,0]
    Vb=UDL_record[n,1]
    Shear=np.zeros(len(X))
    Moment=np.zeros(len(X))
    
    for i,x in enumerate(X):
        
        shear=0
        moment=0
        if x>A:
            shear=shear+Va
            moment=moment-Va*(x-A)
            
        if x>B:
            shear=shear+Vb
            moment=moment-Vb*(x-B)
        
        if x>xstart and x<=xend:
            shear=shear+fy*(x-xstart)
            moment=moment-fy*(x-xstart)*0.5*(x-xstart)
        elif x>xend:
            shear=shear+fy*(xend-xstart)
            moment=moment-fy*(xend-xstart)*(x-xstart-0.5*(xend-xstart))
            
        Shear[i]=shear
        Moment[i]=moment
        
    return Shear , Moment  

    
#Shear and moment calculation due to UVL  
def shear_moment_UVL(n):
    
    xstart=UVL[n,0]
    xend=UVL[n,1]
    fy_start=UVL[n,2]
    fy_end=UVL[n,3]
    Va=UVL_record[n,0]
    Vb=UVL_record[n,1]  
    Shear=np.zeros(len(X))
    Moment=np.zeros(len(X))
    for i,x in enumerate(X):
        shear=0
        moment=0
        if x>A:
            shear=shear+Va
            moment=moment-Va*(x-A)
            
        if x>B:
            shear=shear+Vb
            moment=moment-Vb*(x-B)
        
        if x>xstart and x<=xend:
            if abs(fy_start)>0:
                x_base=x-xstart
                f_cut=fy_start-x_base*(fy_start/(xend-xstart))
                R1=0.5*x_base*(fy_start-f_cut)
                R2=x_base*f_cut
                shear=shear+R1+R2
                moment=moment-R1*(2/3)*x_base-R2*(x_base/2)
            else:
                x_base=x-xstart
                f_cut=fy_end*(x_base/(xend-xstart))
                R=0.5*x_base*f_cut
                shear=shear+R
                moment=moment-R*(x_base/3)
        elif x>xend:
            if abs(fy_start)>0:
                R=0.5*fy_start*(xend-xstart)
                shear=shear+R
                moment=moment-R*(x-(xstart+(1/3)*(xend-xstart)))
            else:
                R=0.5*fy_end*(xend-xstart)
                shear=shear+R
                moment=moment-R*(x-(xstart+(2/3)*(xend-xstart)))
        
        Shear[i]=shear
        Moment[i]=moment        
                       
        
    return Shear , Moment  
        
        
    

     
    


if(nPL>0):
    for n,p in enumerate(pointLoads):
        Shear,Moment=shear_moment_PL(n)
        shearForce=np.append(shearForce,[Shear],axis=0)
        bendingMoment=np.append(bendingMoment,[Moment],axis=0)
        
        
        
if(nPM>0):
    for n,p in enumerate(pointMoments):
        Shear,Moment=shear_moment_PM(n)
        shearForce=np.append(shearForce,[Shear],axis=0)
        bendingMoment=np.append(bendingMoment,[Moment],axis=0)
        
        
if(nUDL>0):
    for n,p in enumerate(UDL):
        Shear,Moment=shear_moment_UDL(n)
        shearForce=np.append(shearForce,[Shear],axis=0)
        bendingMoment=np.append(bendingMoment,[Moment],axis=0)   

if(nUVL>0):
    for n,p in enumerate(UVL):
        Shear,Moment=shear_moment_UVL(n)
        shearForce=np.append(shearForce,[Shear],axis=0)
        bendingMoment=np.append(bendingMoment,[Moment],axis=0)        
        
        
        
        
# Shear Force Diagram
plt.figure(figsize=(10, 4))
plt.plot(X, sum(shearForce), color="green", label="Shear Force")
plt.axhline(0, color='black', linewidth=0.8)
plt.fill_between(X, sum(shearForce), color="green", alpha=0.1)
plt.title("Shear Force Diagram")
plt.xlabel("Distance (m)")
plt.ylabel("Shear Force (kN)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Bending Moment Diagram
plt.figure(figsize=(10, 4))
plt.plot(X, -sum(bendingMoment), color="red", label="Bending Moment")
plt.axhline(0, color='black', linewidth=0.8)
plt.fill_between(X, -sum(bendingMoment), color="red", alpha=0.1)
plt.title("Bending Moment Diagram")
plt.xlabel("Distance (m)")
plt.ylabel("Bending Moment (kNm)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
