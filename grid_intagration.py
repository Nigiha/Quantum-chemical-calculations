import numpy as np

#========================================================================
def grid_integration(func, r_range, pruned, args=(), log_grid=True):
      '''
      func: function to be integrated

      r_range: the radius of grids sphere

      pruned: (r_num, theta_num, phi_num)

      args: parametars

      log_grid: r=b(e^x-1)
      '''

      r_num, theta_num, phi_num = pruned

      r_step = r_range / r_num
      theta_step = np.pi / theta_num
      phi_step = 2 * np.pi / phi_num

      if log_grid==True:

            b=1.0

            x_max=np.log(r_range/b+1.0)
            x_step=x_max/r_num
            x=np.linspace(x_step/2, x_max-x_step/2, r_num)
            r=b*(np.exp(x)-1)
            dr_dx=b*np.exp(x)
            dr=dr_dx*x_step
            dr = dr.reshape(r_num, 1, 1)

            dtheta=theta_step
            dphi=phi_step
      else:
            r=np.linspace(r_step/2, r_range-r_step/2, r_num)
            dr=r_step
            dtheta=theta_step
            dphi=phi_step

      theta=np.linspace(theta_step/2, np.pi-theta_step/2, theta_num)
      phi=np.linspace(phi_step/2, 2*np.pi-phi_step/2, phi_num)

      R, Theta, Phi=np.meshgrid(r, theta, phi, indexing="ij")
        
      dV=R**2*np.sin(Theta)*dr*dtheta*dphi
      func_value=func(R, Theta, Phi, *args)

      result=np.sum(func_value*dV)

      return result

#========================================================================

def func(R, Theta, Phi, a):
      return np.exp(-a*R**2)

a=1.0
k=10
r_range=k/a
pruned=[100, 256, 128]

exact=np.sqrt(np.pi/a)**3
calc=grid_integration(func, r_range, pruned, args=(1.0,), log_grid=True)
error=abs(calc-exact)

print(f"計算結果:{calc:.8f}")
print(f"解析計算:{exact:.8f}")



        
