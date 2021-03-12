// NAME: Zhiguo(Frank) Zhang
// ID: 260550226
package comp559.project;

public class SymplecticEuler implements Integrator {

    @Override
    public String getName() {
        return "symplectic Euler";
    }
    double[] dydt;
    @Override
    public void step(double[] y, int n, double t, double h, double[] yout, Function derivs) {
        // TODO: Objective 7, complete the symplectic Euler integration method.
    	// note you'll need to know how y is packed to properly implement this, so go
    	// look at ParticleSystem.getPhaseSpace()
	    	if ( dydt == null || dydt.length != n ) {
	    		dydt = new double[n];
	    	}
	    	for (int i = 0; i < n; i++)
	    		dydt[i] = 0;
    		derivs.derivs(t, y, dydt);
    		for (int i = 0; i < n; i+=6) {
    			yout[i+3] = y[i+3] + h * dydt[i+3];
    			yout[i+4] = y[i+4] + h * dydt[i+4];
    			yout[i+5] = y[i+5] + h * dydt[i+5];
    			yout[i]   = y[i]   + h * yout[i+3];
    			yout[i+1] = y[i+1] + h * yout[i+4];
    			yout[i+2] = y[i+1] + h * yout[i+5];
    		}

    }

}
