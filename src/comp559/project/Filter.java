// NAME: Zhiguo(Frank) Zhang
// ID: 260550226
package comp559.project;

import no.uib.cipr.matrix.Vector;

/**
 * Velocity filter to use with a conjugate gradients solve
 * @author kry
 */
public interface Filter {

    /**
     * removes disallowed parts of v by projection
     * @param v
     */
    public void filter( Vector v );
    
}
