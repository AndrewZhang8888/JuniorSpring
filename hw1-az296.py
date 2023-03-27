import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.stats as stats
from scipy.stats import pareto

#Question 1
samples=np.linspace(start=50,stop=150,num=1000)

#1.B
#Pareto PDF function
a=3
theta=10

# +
paretos=np.copy(samples)
pareto=lambda x: a*theta**a/(theta+x)**(a+1)
paretos=pareto(paretos)
plt.plot(samples,paretos,label='pareto')
plt.legend()

#1.C
gammas=stats.gamma.pdf(samples,a=1/3,scale=15)
plt.plot(samples,gammas,label='gamma')
plt.legend()
plt.show()

# +
#Question 3
#3.A
gamma_cdf=stats.gamma.cdf(samples,a=1/3,scale=15)
surv=lambda x: 1-x
gamma_surv=surv(gamma_cdf)
gamma_hazard=np.divide(gammas,gamma_surv)
plt.plot(samples,gamma_hazard,label='gamma_hazard')
plt.legend()

#Question 2 in PDF document

pareto_haz=lambda x:a/(theta+x)
pareto_hazard=pareto_haz(samples)
plt.plot(samples,pareto_hazard,label='pareto_hazard')
plt.legend()

#We have found our limits in this question by plotting large x values and their respective hazard functions
large_x=np.full((1),500)
pareto_limit=pareto_haz(large_x)
gamma_limit=stats.gamma.pdf(large_x,a=1/3,scale=15)/stats.gamma.sf(large_x,a=1/3,scale=15)
pareto_limit_line=np.repeat(pareto_limit,1000)
gamma_limit_line=np.repeat(gamma_limit,1000)
plt.plot(samples,gamma_limit_line,label='gamma_limit',linestyle='--',color=(0.12,0.47,0.71))
plt.plot(samples,pareto_limit_line,label='pareto_limit',linestyle='--',color='orange')
plt.legend()
plt.show()
# -

#3.B
'''
gamma_res_life=lambda x: scipy.integrate.quad(stats.gamma.sf(a=1/3,scale=15),x,1000000)
gamma_life=np.array([1])
gamma_life=gamma_res_life(gamma_life)'''
#I unfortunately was limited by scipy in finding an integration of the gamma survival function
#As such, I decided to find the area under the survival curve to the right of d manually and then divide by the sf(d) to find e(d)
def integrate(x,y):
    area=np.trapz(y=y,x=x)
    return area
gamma_life_surv_area=[]
samples1=[]
for i in range(100):
    samples1.append(50+i)
    samples2=np.linspace(start=50+i,stop=1000,num=1000)
    gamma_sf=stats.gamma.sf(samples2,a=1/3,scale=15)
    gamma_single_life=integrate(samples2,gamma_sf)
    gamma_life_surv_area.append(gamma_single_life)
gamma_life_surv=stats.gamma.sf(samples1,a=1/3,scale=15)
gamma_life=np.divide(np.array(gamma_life_surv_area),np.array(gamma_life_surv))
plt.plot(samples1,gamma_life,label='gamma_life')
plt.legend()
plt.show()

#The pareto mean residual life function was provied
pareto_res_life=lambda x:(x+theta)/(a-1)
pareto_life=pareto_res_life(samples)
plt.plot(samples,pareto_life,label='pareto_life')
plt.legend()
plt.show()

# Question 4 in PDF document

# Question 5 in PDF document
