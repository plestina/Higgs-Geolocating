def getGamma33(channel):
    g = 0.038
    if channel==3:  g=0.040
    else: 	    g=0.034 
    print "gamma33(ch={0})={1}".format(channel,g)
    return g

def fa3(k3k1,channel):
    #print "fa3Calculator"
    g=getGamma33(channel)
    return g*k3k1*k3k1/(1+g*k3k1*k3k1)
    
def k3k1(fa3,channel):
    #print "k3k1Calculator"
    #pass
    import math
    g = getGamma33(channel)
    if (1-fa3) > 0 : 
      return math.sqrt(fa3/(g*(1 - fa3) ))
    else :
       raise RuntimeError, "Abs of argument must be < 1 "

       
       
## run the create_RM_cfg() as main()
if __name__ == "__main__":
    values = {}
    values['ch1'] = []
    values['ch2'] = [4.15,12.45]
    values['ch3'] = [5.13,9]
    
    doChannel = 3  #2e2mu
    dofa3=True
    dok3k1=False
    
    for val in values['ch{0}'.format(doChannel)]:
      if dofa3:
	print "k3k1 = {0} ==> f_a3 = {1}".format(val,fa3(val,doChannel))
      elif dok3k1:
	print "f_a3 = {0} ==> k3k1 = {1}".format(val,k3k1(val,doChannel))
	
    
    
