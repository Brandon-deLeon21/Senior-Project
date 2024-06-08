

from mpl_toolkits import mplot3d

from numpy import sin,cos,pi,array,sqrt,arctan,arccos,arcsin,arctan2,zeros
from matplotlib import pyplot as plt

class Planet:
    def __init__(self,name,m,a,e):
        self.name = name
        self.m = m
        self.a = a
        self.e = e
        self.x = [0.]
        self.y = [0.]
        self.z = [0.]
        self.vx = [0.]
        self.vy = [0.]
        self.vz = [0.]
        self.vars = array([0.,0.,0.,0.,0.,0.])
        self.t = [0.]
        #self.kE = [.5*self.m*v**2]
        #self.pE = [- 1.98524*10**(-29)*1.9885*10**30*5.9722*10**24/(self.x[0])]
        self.inOrbit = True
        self.idaOrb = False
        self.earthAp = False
        
    def setupPos(self,x,y,z,vx,vy,vz):
        self.vars = array([x,y,z,vx,vy,vz])
        self.x[0] = self.vars[0]
        self.y[0] = self.vars[1]
        self.z[0] = self.vars[2]
        self.vx[0] = self.vars[3]
        self.vy[0] = self.vars[4]
        self.vz[0] = self.vars[5]
        
    def setupVars(self, r,v,phi,th):
        self.x[0] = r*cos(phi)*sin(pi/2-th)
        self.y[0] = r*sin(phi)*sin(pi/2-th)
        self.z[0] = r*cos(pi/2-th)
        self.vx[0] = v*cos(phi+pi/2)
        self.vy[0] = v*sin(phi+pi/2)
        self.vz[0] = 0
        self.vars = array([self.x[0],self.y[0],self.z[0],self.vx[0],self.vy[0],self.vz[0]])
        print(self.name,self.vars)
        
    def setupV(self,dt,others):
        derv = self.deriv(self.vars,others)
        #print(self.vars)
        #print(derv)
        self.vars[2:] += .5*dt*derv[2:]
        #print(self.vars)
    
    def getApa(self):
        r = sqrt(array(self.x)**2 + array(self.y)**2)
        index = []
        for i in range(len(r)-1):
            if(r[i-1]<r[i] and r[i]>r[i+1]):
                index.append(i)
        apa = []
        apaCom = []
        for i in index:
            apa.append(r[i])
            apaCom.append(array([self.x[i],self.y[i]]))
        return array([index,apa,apaCom])
    
            
    def getPeri(self):
        r = sqrt(array(self.x)**2 + array(self.y)**2)
        index = []#0]
        for i in range(len(r)-1):
            if(r[i-1]>r[i] and r[i]<r[i+1]):
                index.append(i)
        peri = []#r[0]]
        periCom = []#array([self.x[0],self.y[0]])]
        for i in index:
            peri.append(r[i])
            periCom.append(array([self.x[i],self.y[i]]))
        return array([index,peri,periCom])
    
    def findYear(self):
        index = []
        for i in range(len(self.y)-1):
            if(abs(self.y[i-1])>abs(self.y[i])<abs(self.y[i+1]) and self.x[i]>0):
                index.append(i)
                
        return index
        
    def deriv(self,vars,others):
        #isida = False
        x = vars[0]
        y = vars[1]
        z = vars[2]
        vx = vars[3]
        vy = vars[4]
        vz = vars[5]
        Relvx = vx
        Relvy = vy
        Relvz = vz
        G = 1.98524*10**(-29) #AU/yr
        
        dx = vx
        dy = vy
        dz = vz
        
        dvx = 0
        dvy = 0
        dvz = 0

        for body in others:
            r2 = (x-body.x[-1])**2 + (y-body.y[-1])**2 + (z-body.z[-1])**2
            a = -G*body.m/r2  #need to change to G*M
            #print("acceleration",a)
            
            dvx += a*(x-body.x[-1])/sqrt(r2)

            dvy += a*(y-body.y[-1])/sqrt(r2)
            
            dvz += a*(z-body.z[-1])/sqrt(r2)

        #print("dvx",atx,"dvy",aty)
        return array([dx,dy,dz,dvx,dvy,dvz])
    
        
    def leapFrog(self,dt,others):
        """find v*1.5 with old position
        derv = self.deriv(self.vars,others)
        derv[2:] = dt*(.5*derv[2:] + self.vars[2:])"""
        """update position with new velocity"""
        #if at > 0:
            #print(self.vars[2:])
        self.vars[:3] += dt*self.vars[3:]
        k = dt*self.deriv(self.vars,others)
        #print(self.vars)
        #print(k)
        #self.pE.append(- 1.98524*10**(-29)*1.9885*10**30*5.9722*10**24/
                       #sqrt((self.vars[0] - others[0].vars[0])**2 + (self.vars[1] - others[0].vars[1])**2))
        self.vars[3:] += .5*k[3:]
        
        #self.kE.append(.5*self.m*(self.vars[2]**2 + self.vars[3]**2))
        
        self.vars[3:] += .5*k[3:]
        #print(self.vars)
        #print("k:",k)
        """use new pos to find new velocity
        k2 = self.deriv(k,others)
        self.vars[2:] += dt*k2[2:]"""
        self.x.append(self.vars[0])
        self.y.append(self.vars[1])
        self.z.append(self.vars[2])
        self.t.append(self.t[-1]+dt)
        
        self.vx.append(self.vars[3])
        self.vy.append(self.vars[4])
        self.vz.append(self.vars[5])
    
    """def graph(self,s=0,f=-1):
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot3D(self.x[s:f],self.y[s:f],self.z[s:f],',')
        #plt.show()"""
    
    def graphEnergy(self):
        plt.plot(self.t[1:],self.energy[1:])
        plt.show()
        print("here")
        
class Probe(Planet):
    def setupStages(self):
        self.idaApproach = False
        self.phi0 = 0
        self.th0 = 0
        self.stage1 = 1
        self.stage2 = 0
        self.ts1 = 47005
        self.s1dur = 80
        self.ts2 = self.ts1 + self.s1dur+1
        self.s2dur = 195
        self.ts2f = self.ts2+self.s2dur+1
        self.taa = 189716
        self.taadur = 25#30
        self.taa2 = 194689
        self.taadur2 = 335#30
        self.tdropao = 200850
        self.tal = 201480
        self.taldur = 16
        self.tadep = 408250
        self.tadepdur = 667
        self.tea = 509888
        self.teadur = 338
        self.tap = 520003
        self.tapdur = 495
        
    def setupV(self,dt,others,at=0):
        derv = self.deriv(self.vars,others,at)
        #print(self.vars)
        #print(derv)
        self.vars[2:] += .5*dt*derv[2:]
        #print(self.vars)
        
    def deriv(self,vars,others,at,yaw=0,pitch=0):
        #isida = False
        yaw *= pi/180
        pitch *= pi/180
        dphi = 0
        dth = 0
        #if at!=0:
            #print(pitch, yaw)
        x = vars[0]
        y = vars[1]
        z = vars[2]
        vx = vars[3]
        vy = vars[4]
        vz = vars[5]
        Relvx = vx
        Relvy = vy
        Relvz = vz
        G = 1.98524*10**(-29) #AU/yr
        
        dx = vx
        dy = vy
        dz = vz
        
        dvx = 0
        dvy = 0
        dvz = 0

        for body in others:
            r2 = (x-body.x[-1])**2 + (y-body.y[-1])**2 + (z-body.z[-1])**2
            a = -G*body.m/r2  #need to change to G*M
            #print("acceleration",a)
            if body.name == "Earth" and at!=0 and self.inOrbit:# and a>.01:
                Relvx = vx - body.vars[3] - 0.0608575
                Relvy = vy - body.vars[4] - 0.0558099
                Relvz = vz - body.vars[5] + 0.0242669
            if body.name == "Ida" and at!=0 and self.idaOrb:
                Relvx = vx - body.vars[3]
                Relvy = vy - body.vars[4]
                Relvz = vz - body.vars[5]
                
                thpos = arctan(sqrt((x - body.x[-1])**2 + (y - body.y[-1])**2)/(z - body.z[-1]))
                phipos = arctan2(x-body.x[-1],y - body.y[-1])*(-1)+pi/2
                tempv = sqrt(Relvx**2 + Relvy**2 + Relvz**2)
                """print((x-body.x[-1]), + (y-body.y[-1]), (z-body.z[-1]))
                print("a-p V:",sqrt(Relvx**2 + Relvy**2 + Relvz**2))
                print("a-p th:",arctan(sqrt((Relvx)**2 + (Relvy)**2)/(Relvz))*180/pi , thpos*180/pi)
                print("a-p phi:",(arctan2(Relvx,Relvy)*(-1)+pi/2)*180/pi -180 - phipos*180/pi)
                """
            
            dvx += a*(x-body.x[-1])/sqrt(r2)

            dvy += a*(y-body.y[-1])/sqrt(r2)
            
            dvz += a*(z-body.z[-1])/sqrt(r2)
            
        #phi = arctan2(Relvx,Relvy)*(-1)+pi/2
        #th = arctan(sqrt(Relvx**2 + Relvy**2)/Relvz)
        #if at!=0:
            #print("phi",phi*180/pi,"th",th*180/pi)
            
        if self.idaApproach and at!=0:    #pointing self in direction of velocity, need to change idaApp name to velDir or something
                self.th0 = arctan2(sqrt((Relvx)**2 + (Relvy)**2),(Relvz))
                self.phi0 = arctan2(Relvx,Relvy)*(-1)+pi/2
                print("asdfargaerzfgsrthwsdfvs",i)
                #print("ang",self.phi0*180/pi,self.th0*180/pi)
                #print(vx,vy,vz)
                #print(vx/vy)
                
        if self.earthAp:
            print(at)
            print(self.phi0*180/pi,self.th0*180/pi)
        
        if at!=0:
            atx = at*cos(self.phi0 + yaw)*sin(self.th0 + pitch)
            aty = at*sin(self.phi0 + yaw)*sin(self.th0 + pitch)
            atz = at*cos(self.th0 + pitch)
            print(atx,aty,atz)
            
            if self.idaApproach:
                th = arctan2(sqrt((atx)**2 + (aty)**2),(atz))
                phi = arctan2(atx,aty)*(-1)+pi/2
                #print("at",atx,aty,atz)#("ang at",phi*180/pi,th*180/pi)
            
            alpha = 0
            if yaw!=0 or pitch!=0:
                atx1 = at*cos(self.phi0)*sin(self.th0)
                aty1 = at*sin(self.phi0)*sin(self.th0)
                atz1 = at*cos(self.th0)
                
                vxtest = Relvx + atx*dt/2
                vytest = Relvy + aty*dt/2
                vztest = Relvz + atz*dt/2
                
                vxtest1 = Relvx + atx1*dt/2
                vytest1 = Relvy + aty1*dt/2
                vztest1 = Relvz + atz1*dt/2
                
                alpha = arctan2(vxtest,vytest)*(-1)+pi/2
                beta = arctan2(sqrt(vxtest**2 + vytest**2),vztest)
                
                alpha1 = arctan2(vxtest1,vytest1)*(-1)+pi/2
                beta1 = arctan2(sqrt(vxtest1**2 + vytest1**2),vztest1)
                
                dphi = alpha - alpha1
                dth = beta - beta1
                
                #alpha = arccos((vxtest1*vxtest + vytest1*vytest + vztest1*vztest)/(sqrt(vxtest1**2 + vytest1**2 + vztest1**2) * sqrt(vxtest**2 + vytest**2 + vztest**2)))
                #print("DPHI:",dphi*180/pi, dth*180/pi)
                #print("betas",beta1*180/pi,beta*180/pi,alpha*180/pi)
                
            #v0 = sqrt(Relvx**2 + Relvy**2 + Relvz**2)
            #self.th0 += arctan(sqrt((v0*cos(self.phi0)*sin(self.th0) + atx)**2 + (v0*sin(self.phi0)*sin(self.th0) + aty)**2)/(v0*cos(self.th0)+atz))
            #self.phi0 += arctan2((v0*cos(self.phi0)*sin(self.th0) + atx),(v0*sin(self.phi0)*sin(self.th0) + aty))*(-1)+pi/2
            #print(arctan(sqrt((v0*cos(self.phi0)*sin(self.th0) + atx)**2 + (v0*sin(self.phi0)*sin(self.th0) + aty)**2)/(v0*cos(self.th0)+atz)),arctan2((v0*cos(self.phi0)*sin(self.th0) + atx),(v0*sin(self.phi0)*sin(self.th0) + aty))*(-1)+pi/2)
            #print("angles",self.th0*180/pi,self.phi0*180/pi)
            self.phi0 += dphi
            self.th0 += dth
            
            dvx += atx
            dvy += aty
            dvz += atz
        #if at!=0:
            #print(at,sqrt(atx**2 + aty**2 + atz**2))
            #print(atx,aty,atz)
            #print("dvx",atx,"dvy",aty,"dvz",atz,"\n")
            #print(sqrt(Relvx**2 + Relvy**2 + Relvz**2))
        return array([dx,dy,dz,dvx,dvy,dvz])
    
        
    def leapFrog(self,dt,others,at=0,yaw=0,pitch=0):
        """find v*1.5 with old position
        derv = self.deriv(self.vars,others)
        derv[2:] = dt*(.5*derv[2:] + self.vars[2:])"""
        """update position with new velocity"""
        #if at > 0:
            #print(self.vars[2:])
        self.vars[:3] += dt*self.vars[3:]
        k = dt*self.deriv(self.vars,others,at,yaw,pitch)
        #print(self.vars)
        #print(k)
        #self.pE.append(- 1.98524*10**(-29)*1.9885*10**30*5.9722*10**24/
                       #sqrt((self.vars[0] - others[0].vars[0])**2 + (self.vars[1] - others[0].vars[1])**2))
        self.vars[3:] += .5*k[3:]
        
        #self.kE.append(.5*self.m*(self.vars[2]**2 + self.vars[3]**2))
        
        self.vars[3:] += .5*k[3:]
        #print(self.vars)
        #print("k:",k)
        """use new pos to find new velocity
        k2 = self.deriv(k,others)
        self.vars[2:] += dt*k2[2:]"""
        self.x.append(self.vars[0])
        self.y.append(self.vars[1])
        self.z.append(self.vars[2])
        self.t.append(self.t[-1]+dt)
        
        self.vx.append(self.vars[3])
        self.vy.append(self.vars[4])
        self.vz.append(self.vars[5])
        
    def probeDirector(self,i,dt):
        if i<self.ts1:
            self.x.append(earth.vars[0])#+0.00004258750455597226)
            self.y.append(earth.vars[1])
            self.z.append(earth.vars[2])
            self.vx.append(earth.vars[2])
            self.vy.append(earth.vars[3])
            self.vz.append(earth.vars[5])
        if self.ts1<= i<probe.tal + probe.taldur:
            if i == self.ts1:#6243
                self.phi0 = arctan2(earth.vars[0],earth.vars[1])*(-1)+pi/2 + (50.7553+90)*pi/180
                self.th0 = (78.9894)*pi/180
                r = sqrt(earth.vars[0]**2 + earth.vars[1]**2)
                self.vars[0] = r*earth.vars[0]/sqrt(earth.vars[0]**2 + earth.vars[1]**2) + 0.00004258750455597226*cos(self.phi0)*sin(self.th0)
                self.vars[1] = r*earth.vars[1]/sqrt(earth.vars[0]**2 + earth.vars[1]**2) + 0.00004258750455597226*sin(self.phi0)*sin(self.th0)
                self.vars[2] = earth.vars[2] + 0.00004258750455597226*cos(self.th0)
                self.x.append(self.vars[0])
                self.y.append(self.vars[1])
                self.z.append(self.vars[2])
                v = .0000001
                self.vars[3] = v*cos(self.phi0)*sin(self.th0) + earth.vars[3] + 0.0608575#- 3.09705*sin(th0)*10**-9
                self.vars[4] = v*sin(self.phi0)*sin(self.th0) + earth.vars[4] + 0.0558099#+ 3.09705*cos(th0)*10**-9
                self.vars[5] = v*cos(self.th0) + earth.vars[5] - 0.0242669
                print("p init dir:",v*cos(self.phi0),v*sin(self.phi0))
                #print(self.vars[3] - earth.vars[3],self.vars[4] - earth.vars[4])
                probe.setupV(dt,array([earth,sun,jupiter]),76965.3)
                #print(self.vars[3] - earth.vars[3],self.vars[4] - earth.vars[4])
            #start thrusters after 536
            elif self.ts1 < i <(self.ts2):#1670 < i <2400:#was6243-6252
                probe.leapFrog(dt,array([earth,sun,jupiter]),(76965.3+(540410.9-76965.3)/self.s1dur*self.stage1),yaw=0,pitch=1.0)#thrust acceleration equation here already taking into account the change in mass
                self.stage1 += 1
                print(self.stage1-1)
                #print(i, dt, sqrt(self.vars[3]**2 + self.vars[4]**2 + self.vars[5]**2))
            elif (self.ts2) < i <(self.ts2f):
                #probe.inOrbit = False
                probe.leapFrog(dt,array([earth,sun,jupiter]),(25398.2+((121694.4-25398.2)/260)*self.stage2),yaw=7.75,pitch=12.25)#1.35)
                self.stage2 += 1
                #print(stage2)
                #print(i, dt)
            elif i==(self.ts2f):
                probe.leapFrog(dt,array([earth,sun,jupiter]),(25398.2+((121694.4-25398.2)/260)*self.stage2),yaw=7.75,pitch=12.1)#51099.08/2,-15)#47620
                #print(stage2)
                #print(i, dt)
                dt = .00001
                probe.inOrbit = False
                #print((TIME-TimeCount))
                print("launch time:", (TIME - TimeCount)*365.2425*24*3600)
            elif (self.taa) < i <=(self.taa + self.taadur):
                probe.leapFrog(dt,array([earth,sun,jupiter,ida]),18503.7,yaw=-1.5,pitch=2.2)
            elif (self.taa2) < i <(self.taa2 + 93):
                #probe.inOrbit = False
                probe.leapFrog(dt,array([earth,sun,jupiter,ida]),18503.7,yaw=30,pitch=90-self.th0*180/pi)
                #print(i,dt)
            elif (self.taa2+93) < i <(self.taa2 + self.taadur2):
                #probe.inOrbit = False
                probe.leapFrog(dt,array([earth,sun,jupiter,ida]),18503.7,yaw=-14.7,pitch=90-self.th0*180/pi)
                #print(i,dt)
            elif i== (self.taa2 + self.taadur2):
                #probe.inOrbit = False
                probe.leapFrog(dt,array([earth,sun,jupiter,ida]),18503.7,yaw=-14.5,pitch=90-self.th0*180/pi)
                #print(i,dt)
            elif (self.tdropao) < i <(probe.tdropao+2):
                probe.leapFrog(dt,array([earth,sun,jupiter,ida]),-18503.7)
            elif (self.tal) < i <(probe.tal + probe.taldur):
                probe.leapFrog(dt,array([earth,sun,jupiter,ida]),-18503.7,yaw=4,pitch=0)
            elif i<self.taa:
                probe.leapFrog(dt,array([earth,sun,jupiter]))
            else:
                probe.leapFrog(dt,array([earth,sun,jupiter,ida]))
        if probe.tal + probe.taldur<=i<self.tadep:
            self.x.append(ida.vars[0])#+0.00004258750455597226)
            self.y.append(ida.vars[1])
            self.z.append(ida.vars[2])
            self.vx.append(ida.vars[2])
            self.vy.append(ida.vars[3])
            self.vz.append(ida.vars[5])
        if self.tadep <= i:
            if i == self.tadep:
                self.phi0 = arctan2(ida.vars[0],ida.vars[1])*(-1)+pi/2 + (150+90)*pi/180
                #self.th0 = (78.9894)*pi/180
                r = sqrt(ida.vars[0]**2 + ida.vars[1]**2)
                self.vars[0] = r*ida.vars[0]/sqrt(ida.vars[0]**2 + ida.vars[1]**2) + 8.021505*10**-8*cos(self.phi0)*sin(self.th0)
                self.vars[1] = r*ida.vars[1]/sqrt(ida.vars[0]**2 + ida.vars[1]**2) + 8.021505*10**-8*sin(self.phi0)*sin(self.th0)
                self.vars[2] = ida.vars[2] + 8.021505*10**-8*cos(self.th0)
                self.x.append(self.vars[0])
                self.y.append(self.vars[1])
                self.z.append(self.vars[2])
                v = .0000001
                self.vars[3] = v*cos(self.phi0)*sin(self.th0) + ida.vars[3] #+ 0.0608575#- 3.09705*sin(th0)*10**-9
                self.vars[4] = v*sin(self.phi0)*sin(self.th0) + ida.vars[4] #+ 0.0558099#+ 3.09705*cos(th0)*10**-9
                self.vars[5] = v*cos(self.th0) + ida.vars[5] #- 0.0242669
                print("p init dir:",v*cos(self.phi0),v*sin(self.phi0))
                #print(self.vars[3] - earth.vars[3],self.vars[4] - earth.vars[4])
                probe.setupV(dt,array([earth,sun,ida,jupiter]),76965.3)
            elif self.tadep < i < self.tadep+self.tadepdur:
                probe.leapFrog(dt,array([earth,sun,jupiter,ida]),18503.7,yaw=0,pitch=0)
                print(self.phi0*180/pi,self.th0*180/pi)
            elif self.tea < i < self.tea+301:
                probe.leapFrog(dt,array([earth,sun,jupiter,ida]),-18503.7,yaw=-3,pitch=-10)
                print(self.phi0*180/pi,self.th0*180/pi)
            elif self.tea+300 < i < self.tea+self.teadur:
                probe.leapFrog(dt,array([earth,sun,jupiter,ida]),-18503.7,yaw=-1,pitch=-10)
                print(self.phi0*180/pi,self.th0*180/pi)
            elif self.tap < i < self.tap+self.tapdur:
                probe.leapFrog(dt,array([earth,sun,jupiter,ida]),-18503.7,yaw=15,pitch=0)
                print(i)
                #print(self.phi0*180/pi,self.th0*180/pi)
            else:
                probe.leapFrog(dt,array([earth,sun,jupiter,ida]))
        
    
def findInitVelocity(a,e):
    G = 1.98526*10**(-29)
    M = 1.9885*10**30

    return sqrt(G*M*(1+e)/(a*(1-e)))

def graph(title,bodies,s=0,f=-1,tilt=30,rot=-100):
    colors = [',r',',m',',b',',c',',y']
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    i = 0
    for body in bodies:
        #print(body.name)
        ax.plot3D(body.x[s:f],body.y[s:f],body.z[s:f],colors[i])
        i += 1
    ax.set_xlabel('$X(au)$')
    ax.set_ylabel('$Y(au)$')
    ax.set_zlabel('$Z(au)$')
    plt.title(title)
    ax.axes.set_zlim3d(bottom=-.6, top=.6) 
    ax.view_init(tilt, rot)
    plt.show()
    #print("here")
    
def CenteredGraph(title,center,orbital,R,s=0,f=-1,tilt=30,rot=-100):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(array(orbital.x[s:f])-array(center.x[s:f]),array(orbital.y[s:f])-
              array(center.y[s:f]),array(orbital.z[s:f])-array(center.z[s:f]),'.b')
    
    theta = [n*2*pi/100 for n in range(101)] #plots a small circle to represent position of the sun
    xE = array([R*cos(angle) for angle in theta]) 
    yE = array([R*sin(angle) for angle in theta]) 
    z = zeros(101)
    ax.plot3D(xE,yE,z, 'r')
    ax.plot3D(xE,z,yE, 'r')
    ax.plot3D(z,xE,yE, 'r')
    
    ax.set_xlabel('$X$')
    ax.set_ylabel('$Y$')
    ax.set_zlabel('$Z$')
    plt.title(title)
    #ax.axes.set_zlim3d(bottom=-.6, top=.6) 
    ax.view_init(tilt, rot)
    plt.show()
    
    
#def TimeChanges(t):
    
    

"""object atributes"""
#sun initial conditions 
mSol = 1.9885*10**30 #kg
aSol = 3.003*10**-6 #AU
eSol = 0.016983
v0Sol = -1.91873979*10**-5#-1.91875*10**-5#-1.885*10**-5
SIncl = 1.3*pi/180

#jupiter initial conditions
mJupiter = 1.89813*10**27
aJupiter = 5.2029
eJupiter = .0487
JIncl = 1.3*pi/180

#earth initial conditions
mEarth = 5.9722*10**24
aEarth = 1
eEarth = .01671
EIncl = 0

#ida initial conditions
mIda = 3.4*10**16
aIda = 2.861
eIda = 	0.0411
IIncl = 1.132*pi/180

mProbe = 3000
aProbe = 0.0002832326409576363#36000000+6371000
eProbe = 0

"""initializing objects"""
#(mass,semimajor,eccentricity,x0,y0,v0,other object positions)
sun = Planet("Sun",mSol,aSol,eSol) 
jupiter = Planet("Jupiter",mJupiter,aJupiter,eJupiter)
earth = Planet("Earth",mEarth,aEarth,eEarth)
ida = Planet("Ida",mIda,aIda,eIda)
probe = Probe("Probe",mProbe,aProbe,eProbe)

"""(x,y,z,vx,vy,vz)"""
pj = 4.944798694191398
vpj = 2.8907183743443947
phiJ = 294*pi/180
thJ = -1.3*pi/180
#jupiter.setupPos(pj*cos(JIncl),0,pj*sin(JIncl), 0,vpj,0)
jupiter.setupVars(pj,vpj,phiJ,thJ)
"""print(pj*cos(JIncl)+jupiter.vars[0])
print(pj*sin(JIncl)+jupiter.vars[2])
print(vpj+jupiter.vars[4])"""

           
ps = 0.004720075808602222
vps = 0.0027593458727152766
phiS = phiJ-pi
thS = -thJ
#sun.vars = array([ps*cos(SIncl),0,ps*sin(SIncl), 0,vps,0])
sun.setupVars(ps,vps,phiS,thS)


pe = .9880100758
vpe = 6.388897208268812
phiE = 0*pi/180
thE = 0*pi/180
#earth.vars = array([.9880100758,0,sun.vars[2],0,6.388897208268812,sun.vars[5]])
earth.setupVars(pe,vpe,phiE,thE)
print("hwowefn:",earth.vars)
earth.vars[2] = sun.vars[2]
print(earth.vars)
"""print(earth.vars[0]-.9880100758)
print(earth.vars[1])
print(earth.vars[2]- sun.vars[2])
print(earth.vars[3])
print(earth.vars[4]-6.388897208268812)
print(earth.vars[5]-sun.vars[5])"""


pI = 2.74341
vpI = 3.8705450693648475
phiI = 241*pi/180 #110
thI = -1.132*pi/180
#ida.vars = array([-0.9381193488120602+sun.vars[0],2.5774617276146103,-0.05419843127580242,
                  #-3.6364127985602646,-1.3235460181805958+0.0027590699381280052,0.07646587611433431])
ida.setupVars(pI,vpI,phiI,thI)
print("init vel ida:",findInitVelocity(aIda,eIda))

#print(earth.vars[0])
#probe.vars = array([earth.vars[0]+0.0002832326409576363,0,earth.vars[2],0,earth.vars[4]+0.646996,earth.vars[5]])
probe.setupPos(earth.vars[0],earth.vars[1],earth.vars[2],earth.vars[3],earth.vars[4],earth.vars[5])
probe.setupStages()
#print(probe.vars[0])
#probe.setVy(earth.vars[3]+0.646996)
#probe.vx.append(0)
#probe.vy.append(0)


"""solve numerical problem"""
dt = .00001#.0001 #1.753hrs/2
TIME = 0
TimeCount = 0
numBodies = 2
numOrbits = 1#24*6

sun.setupV(dt,array([jupiter]))
jupiter.setupV(dt,array([sun]))
earth.setupV(dt,array([sun,jupiter]))
ida.setupV(dt,array([sun,jupiter,earth]))
probe.setupV(dt,array([earth,sun,jupiter]))
#print(earth.vars[4],probe.vars[4])
#print(probe.vars[0])
ts1 = 10318#6205
ts2 = ts1+9
ts2f = ts2+18#+11
tIdaDeparture = 309880 #this needs to go down to get probe closer to earth

#timeStages = array([ts1,.00001],[ts2,.00001])
for i in range(557100):#10001*numOrbits):#52970#170278
    """Earth launch"""    
    if i == (probe.ts1-1):    #first rocket stage
        dt = .000005
    if i == (probe.ts1):
        dt = .0000001
        TimeCount = TIME
    if i == probe.ts2:
        print("launch time:", (TIME - TimeCount)*365.2425*24*3600)
        TimeCount = TIME
    if i == probe.ts2f:
        dt = .0000001*.2765 #end of second rocket stage
    if i == probe.ts2f+1:
        dt = .00001
    """Asteroid interception"""
    if i == probe.taa-1:
        dt = .0000052
    if i == probe.taa:
        dt = .000001
        probe.idaApproach = True
    if i == probe.taa+probe.taadur:
        dt = .000001*.645
    if i == probe.taa+probe.taadur+1:
        dt = .00001
    if i == 194530:
        dt = .000001
    if i == probe.taa2-1:
        dt = .000001*.82
    if i == probe.taa2:
        dt = .0000001
    if i == (probe.taa2+probe.taadur2):
        dt = .0000001*.6
        probe.idaApproach = False
    if i == (probe.taa2+probe.taadur2+1):
        dt = .0000001
    if i == (probe.tdropao):
        probe.idaOrb = True
        probe.idaApproach = True
    if i == probe.tal-1:
        dt = dt
    if i == probe.tal:
        dt = .00000001
    if i == probe.tal+probe.taldur:
        dt = .00001
    if i == probe.tadep:
        dt = .0000001
    if i == probe.tadep+probe.tadepdur+1:
        dt = .00001
        probe.idaOrb = False
    if i == probe.tea:
        dt = .0000001
    if i == probe.tea+probe.teadur:
        dt = .00001
    if i == probe.tap-1:
        dt = .0000001*.85
    if i == probe.tap:
        dt = .0000001
        probe.earthAp = True
    if i == probe.tap+probe.tapdur:
        dt = .000005
        probe.earthAp = False
    
        
        
    TIME += dt
    sun.leapFrog(dt,array([jupiter]))
    
    jupiter.leapFrog(dt,array([sun]))
    
    earth.leapFrog(dt,array([sun,jupiter]))
    
    ida.leapFrog(dt,array([sun,jupiter,earth]))
    
    #probe.leapFrog(dt,array([earth,sun,jupiter]))
    probe.probeDirector(i,dt)
"""When making the probe part CREATE a function that will intake the start burn times
or durations so that there are not as many if else statements"""
print(probe.vars[0])
m = 0
mi = 0
print(probe.x[0],earth.x[0])
for i in range(len(probe.x)):
    if probe.x[i]-earth.x[i] < m:
        m = probe.x[i]-earth.x[i]
        mi = i
print(m,mi)

"""graph orbits"""
#graph(bodies, startTime, endTime, Tilt, Rotation)
graph(" ",array([sun,jupiter,earth,ida,probe]))
graph(" ",array([sun,jupiter,earth,ida,probe]),tilt=0)
graph(" ",array([sun,jupiter,earth,ida,probe]),tilt=90)

#graph(array([sun,earth,ida]),0,5000,80,-90)

CenteredGraph("Earth Centered",earth,probe,0.0000425875,probe.ts1,probe.ts2f+10)
CenteredGraph("Earth Centered rotated",earth,probe,0.0000425875,probe.ts1,probe.ts2f+10,rot=-170,tilt=5)
CenteredGraph("Earth Centered Top-Down",earth,probe,0.0000425875,probe.ts1,probe.ts2f+10,tilt=90)

s = probe.taa2+probe.taadur2-20
f = 201450
CenteredGraph("Ida Centered",ida,probe,8.021505*10**-8,s,f)
CenteredGraph("Ida Centered",ida,probe,8.021505*10**-8,s,f,tilt=0)
CenteredGraph("Ida Centered",ida,probe,8.021505*10**-8,s,f,tilt=90)
graph("Ast-Pro Side",array([ida,probe]),s,f,tilt=0)
graph("Ast-Pro Top-Down",array([ida,probe]),s,f,tilt=90)

s = probe.tal
f = probe.tal+probe.taldur
CenteredGraph("Ida Landing",ida,probe,8.021505*10**-8,s,f,tilt=90)
graph("probe",array([ida]),s,f,tilt=90)

print("probe v:",sqrt((array(probe.vx[s:f]) - array(ida.vx[s:f]))**2 
                      + array((probe.vy[s:f]) - array(ida.vy[s:f]))**2 
                      + array((probe.vz[s:f]) - array(ida.vz[s:f]))**2))
print("r", (sqrt((array(probe.x[s:f]) - array(ida.x[s:f]))**2 
                      + array((probe.y[s:f]) - array(ida.y[s:f]))**2 
                      + array((probe.z[s:f]) - array(ida.z[s:f]))**2) - 8.021505*10**-8)*149597870700)

print("v prop",(array(probe.vx[s:f]) - array(ida.vx[s:f]))
                      ,array((probe.vy[s:f]) - array(ida.vy[s:f])) 
                      ,array((probe.vz[s:f]) - array(ida.vz[s:f])))

s = probe.tadep
f = probe.tadep+probe.tadepdur
CenteredGraph("Ida departure",ida,probe,8.021505*10**-8,s,f,tilt=90)
CenteredGraph("Ida departure",ida,probe,8.021505*10**-8,s,f,tilt=0)

s = probe.tap+300
f = probe.tap+probe.tapdur+150#520400
colors = [',r',',m',',b',',c',',y']
bodies = array([earth,ida,probe])
ax = plt.axes(projection='3d')
i = 0
for body in bodies:
    ax.plot3D(body.x,body.y,body.z,colors[i])
    i += 1
theta = [n*2*pi/100 for n in range(101)] #plots a small circle to represent position of the sun
xE = array([.08*cos(angle) for angle in theta]) + earth.x[-1]
yE = array([.08*sin(angle) for angle in theta]) + earth.y[-1]
z = zeros(101) + earth.z[-1]
ax.plot3D(xE,yE,z, 'r')
ax.set_xlabel('$X$')
ax.set_ylabel('$Y$')
#plt.title(title)
ax.axes.set_zlim3d(bottom=-.6, top=.6) 
ax.view_init(90, -100)
plt.show()

graph("Earth approach",array([earth,probe]),s,f,tilt=90)
CenteredGraph("Earth Centered",earth,probe,0.0000425875,s,f,tilt=90)
CenteredGraph("Earth Centered",earth,probe,0.0000425875,s,f,tilt=0)
"""
print("probe v:",sqrt((array(probe.vx[s:f]) - array(earth.vx[s:f]))**2 
                      + array((probe.vy[s:f]) - array(earth.vy[s:f]))**2 
                      + array((probe.vz[s:f]) - array(earth.vz[s:f]))**2))"""
"""
graph(array([sun]),tilt=0)
graph(array([ida]),tilt=0)
graph(array([earth]))"""
print("done t:",TIME)#*365.2425)