# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd


#%%

class Physical_parameters:
    def __init__(self):
        self.total_abs=0
        self.total_ref=0
        self.total_trans=0
        self.total_spec=0
    
class Photon_packet:
    def __init__(self, beam_type,r):
        
        if beam_type==1:
            self.x=0;
            self.y=0; 
            
        elif beam_type==2:
            self.x=r*np.random.rand(1,1);
            self.y=0; 
            
        
        elif beam_type==3:
            self.x=r*np.sqrt(-np.log(1-np.random.rand(1,1))/2);
            self.y=0; 
        self.z=0;
        self.ux=0;
        self.uy=0;
        self.uz=1;
        self.w=1; 
        self.t=0; 
        self.re=0; 
        self.out=0; 
        self.max_depth=0;


class Layer:
    def __init__(self,param_array, bottom):
        self.n=param_array[0]
        self.ua=param_array[1]
        self.us=param_array[2]
        self.ut=self.ua+self.us
        self.g=param_array[3]
        self.d=param_array[4]
        self.bottom=bottom
        self.top=self.bottom+self.d
    

#%%

def Scatter(photon,layers):
  
    coszero=(1.0-1.0E-12)
    pi=3.1415926535
    r = 0 + 1.*np.random.rand(1,1)
    if layers.g!=0:
        cos0=((1/(2*layers.g))*(1+layers.g**2-((1-layers.g**2)/(1-layers.g+2*layers.g*r))**2));
        if(cos0<-1):
            cos0=-1
        elif(cos0>1):
            cos0=1;
    else:
        cos0=2*r-1;
          
    sin0=np.sqrt(1-cos0**2)
    r = 0 + 1.*np.random.rand(1,1);
    fi=2*pi*r;
    if photon.uz<-coszero or photon.uz>coszero:
        photon.ux=sin0*np.cos(fi)
        photon.uy=sin0*np.sin(fi)
        photon.uz=np.sign(photon.uz)*cos0
    else:
        temp=np.sqrt(1-photon.uz**2);
        photon.ux=((sin0*(photon.uz*photon.ux*np.cos(fi)-photon.uy*np.sin(fi)))/temp+photon.ux*cos0);
        photon.uy=((sin0*(photon.uy*photon.uz*np.cos(fi)+photon.ux*np.sin(fi)))/temp+photon.uy*cos0);
        photon.uz=(-temp*sin0*np.cos(fi)+photon.uz*cos0);


def MovePhoton(photon,s):
     photon.x=photon.x+photon.ux*s
     photon.y=photon.y+photon.uy*s
     photon.z=photon.z+photon.uz*s

def Calculate_results(physical_parameters,n_photons):
    result=[]
    result.append(physical_parameters.total_spec/n_photons)
    result.append(physical_parameters.total_ref/n_photons)
    result.append(physical_parameters.total_abs/n_photons)
  
    result.append(physical_parameters.total_trans/n_photons)

    return result

def Load_data():
    filepath = input("Enter filepath for layer file: ")
    df=pd.read_csv(filepath)
    layers_info=df.to_numpy()
    layers_list=[]
    curr_height=0
    for i, layer in enumerate(layers_info):
        layers_list.append(Layer(layer,curr_height))
        curr_height=layers_list[i].bottom+layers_list[i].d
    return layers_list

#%%
def roulette(photon):
    if photon.w< 0.0001:
        r = 0 + 1.*np.random.rand(1,1);
        if r<0.1:
            photon.w= photon.w*10;
        else:
            photon.w=0; 

def absorbtion(phy_par,photon,layers):
     if layers.ut!=0:
         delW=(layers.ua/layers.ut)*photon.w   
         photon.w=photon.w-delW
         phy_par.total_abs=phy_par.total_abs+delW


    
    
def enter_skin(photon,layer, phy_par):
    R=((1-layer.n)/(1+layer.n))**2
    photon.w=photon.w-R
    phy_par.total_spec=phy_par.total_spec+R




       

def snelious_law(photon,layers,current_layer,next_layer):
       
      
       ai=np.arccos(abs( photon.uz))  
        
      
       coszero=(1.0-1.0E-12)
       cos90d=1.0E-6;
       tempcos=np.cos(ai);
      
       if layers[current_layer].n>layers[next_layer].n and ai> np.arcsin(layers[next_layer].n/layers[current_layer].n):
           R=1
       elif tempcos>coszero:
           R= ((layers[current_layer].n-layers[next_layer].n)/(layers[current_layer].n+layers[next_layer].n))**2
       elif tempcos<cos90d:
           R=1
       else:
           at= np.arcsin(layers[current_layer].n*np.sin(ai)/layers[next_layer].n)
           if np.sin(at)>=1.0:
               R=1
           else:
               R=1/2*((np.sin(ai-at))**2/(np.sin(ai+at))**2+(np.tan(ai-at))**2/(np.tan(ai+at))**2)
       return R
       
    

def Move_photon_boundry(current_layer,layers,s,d,photon):
    sm=(s-d)*layers[current_layer].ut
    s=d
    photon.x=photon.x+photon.ux*s
    photon.y=photon.y+photon.uy*s
    photon.z=photon.z+photon.uz*s
    return sm,d
    
def mc_loop(layers,n_photons,beam_type,r):
    
    phy_par=Physical_parameters()
    
   
   
    



    n=1
    
   
    n_layers=len(layers)-2

    s=0
  
    j=1
 
    while n<=n_photons:
        
    
        
        photon=Photon_packet(beam_type,r)
        

        hit_bnd=0
        hit_lbnd=0
        hit_hbnd=0
    
        current_layer=1
        next_layer=2

       
        sm=0;
        
        if int(n/n_photons*100)==10*j:
            print(str(10*j)+"% "+"["+(j*4-1)*"="+">"+(40-j*4)*"."+"]")
            j+=1
          
       
        enter_skin(photon,layers[1], phy_par)
   
   
        while photon.w>0 and photon.out==0:
            
            r = 1.0E-10 + 1.*np.random.rand(1,1)
        
            if sm==0:
                sm=-np.log(r)
        
   
            s=sm/layers[current_layer].ut
            sm=0
            
            
            
   
            if photon.uz>0:
                d=(layers[current_layer].top-photon.z)/photon.uz
        
                if s>d:
                    if current_layer==n_layers:
                        hit_hbnd=1
                    else:
                        hit_bnd=1
                    next_layer=current_layer+1
            elif photon.uz<0:
                d=(layers[current_layer].bottom-photon.z)/photon.uz;
                if s>d:
                    if current_layer==1:
                        hit_lbnd=1
                    else:
                        hit_bnd=1
                    next_layer=current_layer-1
   
            
            if hit_hbnd==1:
               sm,d=Move_photon_boundry(current_layer,layers,s,d,photon)
               R=snelious_law(photon, layers,current_layer,next_layer)
               r = 0 + 1.*np.random.rand(1,1);
               if r<=R:
                   photon.uz=-photon.uz
                   
               else:
                   photon.out=1
                   phy_par.total_trans=phy_par.total_trans+photon.w;
               hit_hbnd=0
       
            elif hit_lbnd==1:
              
               sm,d=Move_photon_boundry(current_layer,layers,s,d,photon)
               R=snelious_law(photon, layers,current_layer,next_layer)
               r = 0 + 1.*np.random.rand(1,1);
               if r<=R:
                   photon.uz=-photon.uz
             
               else:
                   photon.out=1
                   phy_par.total_ref=phy_par.total_ref+photon.w;
               hit_lbnd=0;
            elif hit_bnd==1:
               sm,d=Move_photon_boundry(current_layer,layers,s,d,photon)
               R=snelious_law(photon, layers,current_layer,next_layer)
               r = 0 + 1.*np.random.rand(1,1);
               if r<=R:
                   photon.uz=-photon.uz
                 
               else:
                   ai=np.arccos(abs( photon.uz))
                   at= np.arcsin(layers[current_layer].n*np.sin(ai)/layers[next_layer].n)
                   photon.ux=  photon.ux*layers[current_layer].n/layers[next_layer].n
                   photon.uy=photon.uy*layers[current_layer].n/layers[next_layer].n
                   photon.uz=np.sign(photon.uz)*np.cos(at);
            
           
                   current_layer=next_layer
                         
               hit_bnd=0;
           
    
            else:
   
 
              
              
          
               MovePhoton(photon,s)
               absorbtion(phy_par, photon, layers[current_layer])
         

               Scatter(photon,layers[current_layer])

   
               
               roulette(photon)
            if photon.w==0 or photon.out==1:
                   n=n+1
    results=Calculate_results(phy_par,n_photons)
    return results





def main():
    results=[]
    r=0
    try:
        layers=Load_data()
    except:
        print("invalid filepath")
        return
    try:
        beam_type=int(input("Enter type of beam: \n 1 - Pencil \n 2 - Flat-field \n 3 - Gaussian \n"))
        assert beam_type==1 or beam_type==2 or beam_type==3
    except:
        print("Unknown type")
        return
    if beam_type==2 or beam_type==3:
        try:
            r=int(input("Enter radius of beam:"))
            assert r>0
        except:
            print("Radius must be a positive integer")
    
    try:
        n_photons=int(input("How many photon packets: "))
        assert n_photons>0
    except:
        print("Photons number must be a positive integer")
        return
    try:
        n_sim=int(input("How many runs: "))
        assert n_sim>0
    except:
        print("Simulations number must be a positive integer")
        return
    for i in range(n_sim):
        results.append(mc_loop(layers, n_photons,beam_type,r))
    print("1.Specular reflectance \n 2.Diffuse reflectance \n 3. Absorbtion \n 4. Transmitance \n "+results)
if __name__ == "__main__": main()