// NAME: Zhiguo(Frank) Zhang
// ID: 260550226
package comp559.project;

import java.util.LinkedList;
import java.util.List;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;

import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.VerticalFlowPanel;
import mintools.viewer.EasyViewer;
import mintools.viewer.SceneGraphNode;

/**
 * Implementation of a simple particle system
 * @author kry
 */
public class ParticleSystem implements SceneGraphNode, Function, Filter {

	public List<Particle> particles = new LinkedList<Particle>();
	public List<Particle> objParticles = new LinkedList<Particle>();
	
	public List<Spring> springs = new LinkedList<Spring>();
	public List<Spring> objSprings = new LinkedList<Spring>();

	public List<Triangle> triangles = new LinkedList<Triangle>();
	public List<Triangle> objTriangles = new LinkedList<Triangle>();

	private Particle[][] clothMesh = null;

	private Vector3d[][] meshNormals = null;
	
	private int index;
	private int recursionLevel = 2;
	
//	private RobustCCD robustCCD = new RobustCCD();
	/**
	 * Creates an empty particle system
	 */
	public ParticleSystem() {
		createSystem();
	}

	/**
	 * Creates one of a number of simple test systems.
	 * @param which
	 */
	public void createSystem() {
		clearParticles();
		createCloth();
		updateObjLists();
		init();
	}

	private void updateObjLists() {
		this.index = 0;
		double t = (1.0 + Math.sqrt(5.0)) / 2.0;
		addVertex(new Point3d(-1, t, 0));
        addVertex(new Point3d(1, t, 0));
        addVertex(new Point3d(-1, -t, 0));
        addVertex(new Point3d(1, -t, 0));

        addVertex(new Point3d(0, -1, t));
        addVertex(new Point3d(0, 1, t));
        addVertex(new Point3d(0, -1, -t));
        addVertex(new Point3d(0, 1, -t));

        addVertex(new Point3d(t, 0, -1));
        addVertex(new Point3d(t, 0, 1));
        addVertex(new Point3d(-t, 0, -1));
        addVertex(new Point3d(-t, 0, 1));
		
		addTriangle(0, 11, 5);
		addTriangle(0, 5, 1);
		addTriangle(0, 1, 7);
		addTriangle(0, 7, 10);
		addTriangle(0, 10, 11);
		
		addTriangle(1, 5, 9);
		addTriangle(5, 11, 4);
		addTriangle(11, 10, 2);
		addTriangle(10, 7, 6);
		addTriangle(7, 1, 8);

		addTriangle(3, 9, 4);
		addTriangle(3, 4, 2);
		addTriangle(3, 2, 6);
		addTriangle(3, 6, 8);
		addTriangle(3, 8, 9);

		addTriangle(4, 9, 5);
		addTriangle(2, 4, 11);
		addTriangle(6, 2, 10);
		addTriangle(8, 6, 7);
		addTriangle(9, 8, 1);
		recursionLevel = sphereSubdivision.getValue();
		for (int i = 0; i < recursionLevel; i++)
        {
            List<Triangle> objTriangles2 = new LinkedList<Triangle>();
            for (Triangle tri : objTriangles)
            {
                // replace triangle by 4 triangles
                int a = getMiddlePoint(tri.p1.index, tri.p2.index);
                int b = getMiddlePoint(tri.p2.index, tri.p3.index);
                int c = getMiddlePoint(tri.p3.index, tri.p1.index);

                objTriangles2.add(new Triangle(objParticles.get(tri.p1.index), objParticles.get(a), objParticles.get(c)));
                objTriangles2.add(new Triangle(objParticles.get(tri.p2.index), objParticles.get(b), objParticles.get(a)));
                objTriangles2.add(new Triangle(objParticles.get(tri.p3.index), objParticles.get(c), objParticles.get(b)));
                objTriangles2.add(new Triangle(objParticles.get(a), objParticles.get(b), objParticles.get(c)));
            }
            objTriangles = objTriangles2;
        }
		for (Triangle tri : objTriangles) {
			addSpring(tri.p1.index, tri.p2.index);
			addSpring(tri.p1.index, tri.p3.index);
			addSpring(tri.p2.index, tri.p3.index);
		}
		List<Spring> removeList = new LinkedList<Spring>();
		for (Spring s1 : objSprings) {
			for (Spring s2 : objSprings) {
				if (s1 != s2 && ((s1.p1 == s2.p1 && s1.p2 == s2.p2) || (s1.p1 == s2.p2 && s1.p2 == s2.p1))) {
					removeList.add(s2);
				}
			}
		}
		for (Spring s : removeList) {
			s.p1.springs.remove(s);
			s.p2.springs.remove(s);
			objSprings.remove(s);
		}
	}
	
	private int addVertex(Point3d pt) {
		Particle p = new Particle(pt.x,pt.y,pt.z, 0,0,0,objParticles.size());
		p.p.scale(0.9/p.p.distance(new Point3d(0,0,0)));
		p.pinned = true;
		objParticles.add(p);
		return index++;
	}
	
	private void addSpring(int i, int j) {
		objSprings.add(new Spring(objParticles.get(i), objParticles.get(j)));
	}
	
	private void addTriangle(int i, int j, int k) {
		objTriangles.add(new Triangle(objParticles.get(i), objParticles.get(j), objParticles.get(k)));
	}
	
	private int getMiddlePoint(int p1, int p2) {
        Point3d point1 = objParticles.get(p1).p;
        Point3d point2 = objParticles.get(p2).p;
        Point3d middle = new Point3d();
        middle.add(point1, point2);
        middle.scale(0.5);
        return addVertex(middle);
	}
	/**
	 * Create a rectangular cloth
	 */
	private void createCloth() {
		int dim = (int) complexity.getValue();
		clothMesh = new Particle[dim][dim];        
		meshNormals = new Vector3d[dim][dim];

		// creating vertices
		for ( int i = 0; i < dim; i++ ) {
			for ( int j = 0; j < dim; j++ ) {               
				clothMesh[i][j] = new Particle(-2 + 4.0 * i/(dim-1), 2.5 - 1.0 * j /(dim-1), -2 + 4.0 * j/(dim-1), 0,0,0, particles.size() );
				particles.add(clothMesh[i][j]);                
				clothMesh[i][j].mass = 1.0 / (dim * dim);                    
				meshNormals[i][j] = new Vector3d();
			}
		}

		// connecting each particle with 8 neighbors
		for ( int i = 0; i < dim; i++ ) {
			for ( int j = 0; j < dim; j++ ) {
				if ( i < dim - 1 ) {
					if (j == 0 || j == dim - 1) springs.add( new Spring( clothMesh[i][j], clothMesh[i+1][j], true) );
					else springs.add( new Spring( clothMesh[i][j], clothMesh[i+1][j]) );
				}
				if ( j < dim - 1 ) {
					if (i == 0 || i == dim - 1) springs.add( new Spring( clothMesh[i][j], clothMesh[i][j+1], true) );
					else springs.add( new Spring( clothMesh[i][j], clothMesh[i][j+1]) );
				}
				if ( i < dim - 1 && j < dim - 1 ) {
					springs.add( new Spring( clothMesh[i][j], clothMesh[i+1][j+1]) );
					springs.add( new Spring( clothMesh[i+1][j], clothMesh[i][j+1]) );
				}
			}                
		}

		// add all triangles
		for ( int i = 0; i < dim - 1; i++ ) {
			for ( int j = 0; j < dim - 1; j++ ) {
				triangles.add(new Triangle(clothMesh[i][j], clothMesh[i+1][j], clothMesh[i+1][j+1]));
				triangles.add(new Triangle(clothMesh[i][j], clothMesh[i+1][j+1], clothMesh[i][j+1]));
			}
		}

		// pinned at two corners
		if (pinTwoPoint.getValue()) {
			clothMesh[0][0].pinned = true;
			clothMesh[dim-1][0].pinned = true;
		} else if (pinTop.getValue()) {
			for(int i = 0; i < dim; i++) {
				clothMesh[i][0].pinned = true;
			}
		}
	}
	/**
	 * Gets the particles in the system
	 * @return the particle set
	 */
	public List<Particle> getParticles() {
		return particles;
	}

	/**
	 * Gets the springs in the system
	 * @return the spring list
	 */
	public List<Spring> getSprings() {
		return springs;
	}

	/**
	 * Resets the positions of all particles
	 */
	public void resetParticles() {
		for ( Particle p : particles ) {
			p.reset();
		}
		time = 0;
	}

	/**
	 * Deletes all particles
	 */
	public void clearParticles() {
		particles.clear();
		objParticles.clear();
		springs.clear();
		objSprings.clear();
		triangles.clear();
		objTriangles.clear();
		clothMesh = null;
		meshNormals = null;
	}

	/**
	 * Gets the phase space state of the particle system
	 * @param y
	 */
	public void getPhaseSpace( double[] y ) {
		int count = 0;
		for ( Particle p : particles ) {
			y[count++] = p.p.x;
			y[count++] = p.p.y;
			y[count++] = p.p.z;
			y[count++] = p.v.x;
			y[count++] = p.v.y;
			y[count++] = p.v.z;
		}
	}

	/**
	 * Gets the dimension of the phase space state
	 * (particles * 3 dimensions * 2 for velocity and position)
	 * @return dimension
	 */
	public int getPhaseSpaceDim() {        
		return particles.size() * 6;
	}

	/**
	 * Sets the phase space state of the particle system
	 * @param y
	 */
	public void setPhaseSpace( double[] y ) {
		int count = 0;
		for ( Particle p : particles ) {
			if ( p.pinned ) {
				count += 6;
			} else {
				p.p.x = y[count++];
				p.p.y = y[count++];
				p.p.z = y[count++];
				p.v.x = y[count++];
				p.v.y = y[count++];
				p.v.z = y[count++];
			}
		}
	}


	/** Elapsed simulation time */
	public double time = 0;

	/** The explicit integrator to use, if not performing backward Euler implicit integration */
	public Integrator integrator;

	public double[] state = new double[1];
	public double[] stateOut = new double[1];

	// these get created in init() and are probably useful for Backward Euler computations
	private ConjugateGradientMTJ CG;
	private DenseMatrix A;
	private DenseMatrix dfdx;
	private DenseMatrix dfdv;
	private DenseVector deltaxdot;
	private DenseVector b;
	private DenseVector f;
	private DenseVector xdot;
	private DenseVector x0, v0;

	/**
	 * Initializes the system 
	 * Allocates the arrays and vectors necessary for the solve of the full system
	 */
	public void init() {
		int N = particles.size();
		// create matrix and vectors for solve
		CG = new ConjugateGradientMTJ(3*N);
		CG.setFilter(this);
		A = new DenseMatrix(3*N, 3*N);
		dfdx = new DenseMatrix(3*N, 3*N);
		dfdv = new DenseMatrix(3*N, 3*N);
		deltaxdot = new DenseVector(3*N);
		b = new DenseVector(3*N);
		f = new DenseVector(3*N);
		xdot = new DenseVector(3*N);
		x0 = new DenseVector(3*N);
		v0 = new DenseVector(3*N);
	}

	/**
	 * Fills in the provided vector with the particle velocities.
	 * @param xd
	 */
	//	private void getVelocities(DenseVector xd) {
	//		for ( Particle p : particles ) {
	//			int j = p.index * 3;
	//			if( p.pinned ) {
	//				xd.set( j,   0 );
	//				xd.set( j+1, 0 );
	//				xd.set( j+2, 0 );
	//			} else {
	//				xd.set( j,   p.v.x );
	//				xd.set( j+1, p.v.y );
	//				xd.set( j+2, p.v.z );
	//			}
	//		}       
	//	}

	/**
	 * Sets the velocities of the particles given a vector
	 * @param xd
	 */
	//	private void setVelocities(DenseVector xd) {
	//		for ( Particle p : particles ) {
	//			int j = p.index * 3;
	//			if( p.pinned ) {
	//				p.v.set(0,0,0);
	//			} else {
	//				p.v.x = xd.get(j);
	//				p.v.y = xd.get(j+1);
	//				p.v.z = xd.get(j+2);
	//			}
	//		}
	//	}

	/**
	 *  Evaluates derivatives for ODE integration.
	 * @param t time 
	 * @param y state
	 * @param dydt to be filled with the derivative
	 */
	@Override
	public void derivs(double t, double[] y, double[] dydt) {
		// set particle positions to given values
		setPhaseSpace( y );

		// TODO: Objective 2, for explicit integrators, compute forces, and accelerations, and set dydt
		for (Particle p : particles) {
			p.clearForce();
			if (useGravity.getValue()) {
				p.addForce(new Vector3d(0, - p.mass * gravity.getValue(), 0));
			}
			p.addForce(new Vector3d( - viscousDamping.getValue() * p.v.x,  - viscousDamping.getValue() * p.v.y,  - viscousDamping.getValue() * p.v.z));
		}
		for (Spring s : springs) {
			s.apply();
		}
		int count = 0;
		for ( Particle p : particles ) {
			dydt[count++] = p.v.x;
			dydt[count++] = p.v.y;
			dydt[count++] = p.v.z;
			dydt[count++] = p.f.x / p.mass;
			dydt[count++] = p.f.y / p.mass;
			dydt[count++] = p.f.z / p.mass;
		}
	}

	/** Time in seconds that was necessary to advance the system */
	public double computeTime;

	/**
	 * Advances the state of the system
	 * @param elapsed
	 */
	public boolean advanceTime( double elapsed ) {
		boolean resultOK = true;
		Spring.k = springStiffness.getValue();
		Spring.c = springDamping.getValue();

		int n = getPhaseSpaceDim();

		long now = System.nanoTime();      
		
//		robustCCD.restitution = restitution.getValue();
//        robustCCD.MAX_ITERATION = iterations.getValue();
//        robustCCD.H = H.getValue();
        
		if ( explicit.getValue() ) {
			if ( n != state.length ) {
				state = new double[n];
				stateOut = new double[n];
			}
			getPhaseSpace(state);         
			integrator.step( state, n, time, elapsed, stateOut, this);
			if ( collision.getValue() ) {
//	            if ( ! robustCCD.check( elapsed, this ) ) {
//	                resultOK = false;
//	            }
				iterativeCheck(elapsed, 0);
	        }
			setPhaseSpace(stateOut);
		} else {        
			if ( f == null || f.size() != n ) {
				init();
			}

			// TODO: Objective 8, your backward Euler implementation will go here!
			// Note that the init() method called above creates a bunch of very 
			// useful MTJ working variables for you, and the ConjugateGradientMTJ object.
			// Go look at that code now!
			if ( n != state.length ) {
				state = new double[n];
				stateOut = new double[n];
			}
			getPhaseSpace(state);         
			implicitEuler(state, elapsed, stateOut);
			if ( collision.getValue() ) {
//	            if ( ! robustCCD.check( elapsed, this ) ) {
//	                resultOK = false;
//	            }
				iterativeCheck(elapsed, 0);
	        }
			setPhaseSpace(stateOut);            
		}
		time = time + elapsed;

		if(repulsion.getValue()) repulsion(elapsed);

		if(usePostStep.getValue()) {
			if(useSphere1.getValue()) postStep(1);
			if(useSphere2.getValue()) postStep(2);
		}
		computeTime = (System.nanoTime() - now) / 1e9;
		return resultOK;
	}

	public void implicitEuler(double[] y, double h, double[] yout) {
		for (int i = 0; i < y.length; i+=6) {
			x0.set(i/2, y[i]);
			x0.set(i/2 + 1, y[i+1]);
			x0.set(i/2 + 2, y[i+2]);
			v0.set(i/2, y[i+3]);
			v0.set(i/2 + 1, y[i+4]);
			v0.set(i/2 + 2, y[i+5]);
		}

		for (Particle p : particles) {
			f.add(3 * p.index, -viscousDamping.getValue() * p.v.x);
			f.add(3 * p.index + 1, -viscousDamping.getValue() * p.v.y);
			f.add(3 * p.index + 2, -viscousDamping.getValue() * p.v.z);
			if (useGravity.getValue()) f.add(3 * p.index + 1, - p.mass * gravity.getValue());
		}

		for (Spring s : springs) {
			s.addForce(f);
			s.addDfdx(dfdx);
			s.addDfdv(dfdv);
		}

		// Compute b
		b.set(dfdx.multAdd(h, v0, f));
		b.scale(h);

		// Compute A
		A.zero();
		for (Particle p : particles) {
			A.set(3 * p.index,     3 * p.index,     1.0); // maybe set to mass 
			A.set(3 * p.index + 1, 3 * p.index + 1, 1.0); // maybe set to mass
			A.set(3 * p.index + 2, 3 * p.index + 2, 1.0); // maybe set to mass
		}
		dfdv.scale(-h);
		dfdx.scale(-h * h);
		A.add(dfdv);
		A.add(dfdx);

		// Solve
		CG.solve(A, b, deltaxdot, 10*iterations.getValue());

		xdot.set(v0.add(deltaxdot));
		x0.add(xdot.scale(h));
		for (int i = 0; i < yout.length; i+=6) {
			yout[i]   = x0.get(i/2);
			yout[i+1] = x0.get(i/2 + 1);
			yout[i+2] = x0.get(i/2 + 2);
			yout[i+3] = v0.get(i/2);
			yout[i+4] = v0.get(i/2 + 1);
			yout[i+5] = v0.get(i/2 + 2);
		}
	}

	@Override
	public void filter(Vector v) {
		for ( Particle p : particles ) {
			if ( !p.pinned ) continue;
			v.set( p.index*3+0, 0 );
			v.set( p.index*3+1, 0 );
			v.set( p.index*3+2, 0 );
		}
	}

	/**
	 * Creates a new particle and adds it to the system
	 * @param x
	 * @param y
	 * @param vx
	 * @param vy
	 * @return the new particle
	 */
	public Particle createParticle( double x, double y, double z, double vx, double vy, double vz ) {
		Particle p = new Particle( x, y, z, vx, vy, vz, particles.size() );
		particles.add( p );
		return p;
	}

	public void remove( Particle p ) {
		for ( Spring s : p.springs ) {
			Particle other = s.p1 == p ? s.p2 : s.p1; 
			other.springs.remove( s );
			springs.remove( s );
		}
		p.springs.clear(); // not really necessary
		particles.remove( p );
		// reset indices of each particle :(
		for ( int i = 0 ; i < particles.size(); i++ ) {
			particles.get(i).index = i;
		}
	}

	/**
	 * Creates a new spring between two particles and adds it to the system.
	 * @param p1
	 * @param p2
	 * @return the new spring
	 */
	public Spring createSpring( Particle p1, Particle p2 ) {
		Spring s = new Spring( p1, p2 ); 
		springs.add( s );         
		return s;
	}

	/**
	 * Removes a spring between p1 and p2 if it exists, does nothing otherwise
	 * @param p1
	 * @param p2
	 * @return true if the spring was found and removed
	 */
	public boolean removeSpring( Particle p1, Particle p2 ) {
		Spring found = null;
		for ( Spring s : springs ) {
			if ( ( s.p1 == p1 && s.p2 == p2 ) || ( s.p1 == p2 && s.p2 == p1 ) ) {
				found = s;
				break;
			}
		}
		if ( found != null ) {
			found.p1.springs.remove(found);
			found.p2.springs.remove(found);
			springs.remove(found);
			return true;
		}
		return false;
	}

	@Override
	public void init(GLAutoDrawable drawable) {
		// do nothing
	}

	@Override
	public void display(GLAutoDrawable drawable) {
		GL2 gl = drawable.getGL().getGL2();
		gl.glDisable(GL2.GL_LIGHTING);
		if (drawVertices.getValue()) {
			gl.glPointSize( 5 );
			gl.glBegin( GL.GL_POINTS );
			for ( Particle p : particles ) {
				double alpha = 0.5;
				if ( p.pinned ) {
					gl.glColor4d( 1, 0, 0, alpha );
				} else {
					gl.glColor4d( p.color.x, p.color.y, p.color.z, alpha );
				}
				gl.glVertex3d( p.p.x, p.p.y, p.p.z);
			}
			for ( Particle p : objParticles ) {
				double alpha = 0.5;
				if ( p.pinned ) {
					gl.glColor4d( 1, 0, 0, alpha );
				} else {
					gl.glColor4d( p.color.x, p.color.y, p.color.z, alpha );
				}
				gl.glVertex3d( p.p.x, p.p.y, p.p.z);
			}
			gl.glEnd();
		}

		gl.glEnable(GL2.GL_LIGHTING);
		displayCloth(drawable);
		displayObj(drawable);
		
		if ( useSphere1.getValue() ) {
			float[] frontColour = { .1f, .4f, .9f, 1};
			float[] specColour  = { .9f, .4f, .1f, 1};
			gl.glMaterialfv( GL.GL_FRONT, GL2.GL_DIFFUSE, frontColour, 0 );
			gl.glMaterialfv( GL.GL_FRONT, GL2.GL_SPECULAR, specColour, 0 );
			gl.glMateriali( GL.GL_FRONT, GL2.GL_SHININESS, 92 );
			gl.glPushMatrix();
			gl.glTranslated( x1pos.getValue(), y1pos.getValue(), z1pos.getValue() );
			EasyViewer.glut.glutSolidSphere( radius1.getValue()*radiusratio.getValue(), 128, 64 );
			gl.glPopMatrix();
		}

		gl.glDisable(GL2.GL_LIGHTING);

		if ( drawGrid.getValue() ) {
			gl.glColor4d(0,.5,.5,.5);
			gl.glLineWidth(2f);
			gl.glBegin( GL.GL_LINES );
			for (Spring s : springs) {
				gl.glVertex3d( s.p1.p.x, s.p1.p.y, s.p1.p.z );
				gl.glVertex3d( s.p2.p.x, s.p2.p.y, s.p2.p.z );
			}
			gl.glEnd();
		}
	}

	/**
	 * display cloth with interpolated shading or flat shading
	 * @param drawable
	 */
	private void displayCloth(GLAutoDrawable drawable) {
		GL2 gl = drawable.getGL().getGL2();
		if ( clothMesh != null && drawCloth.getValue() ) {
			Vector3d v1 = new Vector3d();
			Vector3d v2 = new Vector3d();
			gl.glDisable( GL.GL_CULL_FACE );
			gl.glLightModeli( GL2.GL_LIGHT_MODEL_TWO_SIDE, GL.GL_TRUE );            
			gl.glMaterialfv( GL.GL_FRONT, GL2.GL_DIFFUSE, new float[]{   0,  .9f,  .3f, 1}, 0 );
			gl.glMaterialfv( GL.GL_BACK,  GL2.GL_DIFFUSE, new float[]{ .9f,    0,    0, 1}, 0 );

			if ( !flatShaded.getValue() ) {
				for ( int i = 0; i < clothMesh.length; i++ ) {
					for ( int j = 0; j < clothMesh[0].length; j++ ) {
						if ( i > 0 && i < clothMesh.length-1) {
							v1.sub( clothMesh[i+1][j].p, clothMesh[i-1][j].p );                        
							v1.scale(0.5);
						} else if ( i > 0 ) v1.sub( clothMesh[i][j].p, clothMesh[i-1][j].p );                        
						else v1.sub( clothMesh[i+1][j].p, clothMesh[i][j].p );
						if ( j > 0 && j < clothMesh[0].length-1) {
							v2.sub( clothMesh[i][j+1].p, clothMesh[i][j-1].p );                        
							v2.scale(0.5);
						} else if ( j > 0 ) v2.sub( clothMesh[i][j].p, clothMesh[i][j-1].p );                        
						else v2.sub( clothMesh[i][j+1].p, clothMesh[i][j].p );
						meshNormals[i][j].cross(v1, v2);
						meshNormals[i][j].normalize();
					}
				}
				gl.glBegin(GL.GL_TRIANGLES);
				for ( int i = 0; i < clothMesh.length-1; i++ ) {
					for ( int j = 0; j < clothMesh[0].length-1; j++ ) {
						Point3d p1 = clothMesh[i][j].p;
						Point3d p2 = clothMesh[i+1][j].p;
						Point3d p3 = clothMesh[i+1][j+1].p;
						gl.glNormal3d( meshNormals[i][j].x,meshNormals[i][j].y,meshNormals[i][j].z);
						gl.glVertex3d(p1.x,p1.y,p1.z);
						gl.glNormal3d( meshNormals[i+1][j].x,meshNormals[i+1][j].y,meshNormals[i+1][j].z);
						gl.glVertex3d(p2.x,p2.y,p2.z);
						gl.glNormal3d( meshNormals[i+1][j+1].x,meshNormals[i+1][j+1].y,meshNormals[i+1][j+1].z);
						gl.glVertex3d(p3.x,p3.y,p3.z);                  
						p1 = clothMesh[i][j].p;
						p2 = clothMesh[i+1][j+1].p;
						p3 = clothMesh[i][j+1].p;                            
						gl.glNormal3d( meshNormals[i][j].x,meshNormals[i][j].y,meshNormals[i][j].z);
						gl.glVertex3d(p1.x,p1.y,p1.z);
						gl.glNormal3d( meshNormals[i+1][j+1].x,meshNormals[i+1][j+1].y,meshNormals[i+1][j+1].z);
						gl.glVertex3d(p2.x,p2.y,p2.z);
						gl.glNormal3d( meshNormals[i][j+1].x,meshNormals[i][j+1].y,meshNormals[i][j+1].z);
						gl.glVertex3d(p3.x,p3.y,p3.z);                    
					}
				}
				gl.glEnd();
			} else {                
				Vector3d n = new Vector3d();                
				gl.glBegin(GL.GL_TRIANGLES);
				for ( int i = 0; i < clothMesh.length-1; i++ ) {
					for ( int j = 0; j < clothMesh[0].length-1; j++ ) {
						Point3d p1 = clothMesh[i][j].p;
						Point3d p2 = clothMesh[i+1][j].p;
						Point3d p3 = clothMesh[i+1][j+1].p;
						v1.sub(p2,p1);
						v2.sub(p3,p1);
						n.cross(v1,v2);
						n.normalize();
						gl.glNormal3d(n.x,n.y,n.z);
						gl.glVertex3d(p1.x,p1.y,p1.z);
						gl.glVertex3d(p2.x,p2.y,p2.z);
						gl.glVertex3d(p3.x,p3.y,p3.z);
						p1 = clothMesh[i][j].p;
						p2 = clothMesh[i+1][j+1].p;
						p3 = clothMesh[i][j+1].p;
						v1.sub(p2,p1);
						v2.sub(p3,p1);
						n.cross(v1,v2);
						n.normalize();
						gl.glNormal3d(n.x,n.y,n.z);
						gl.glVertex3d(p1.x,p1.y,p1.z);
						gl.glVertex3d(p2.x,p2.y,p2.z);
						gl.glVertex3d(p3.x,p3.y,p3.z);
					}
				}
				gl.glEnd();
			}            
			gl.glLightModeli( GL2.GL_LIGHT_MODEL_TWO_SIDE, GL.GL_FALSE );
			gl.glEnable( GL.GL_CULL_FACE );
		}
	}
	
	
	private void displayObj(GLAutoDrawable drawable) {
		GL2 gl = drawable.getGL().getGL2();
		if ( /*objMesh != null &&*/ drawObj.getValue() ) {
			Vector3d v1 = new Vector3d();
			Vector3d v2 = new Vector3d();
			gl.glDisable( GL.GL_CULL_FACE );
			gl.glLightModeli( GL2.GL_LIGHT_MODEL_TWO_SIDE, GL.GL_TRUE );            
			gl.glMaterialfv( GL.GL_FRONT, GL2.GL_DIFFUSE, new float[]{ .9f,  .9f,    0, 1}, 0 );
			gl.glMaterialfv( GL.GL_BACK,  GL2.GL_DIFFUSE, new float[]{ .9f,    0,  .9f, 1}, 0 );

			               
			Vector3d n = new Vector3d();                
			gl.glBegin(GL.GL_TRIANGLES);
			for (Triangle t : objTriangles) {
				Point3d p1 = t.p1.p;
				Point3d p2 = t.p2.p;
				Point3d p3 = t.p3.p;
				v1.sub(p2,p1);
				v2.sub(p3,p1);
				n.cross(v1,v2);
				n.normalize();
				gl.glNormal3d(n.x,n.y,n.z);
				gl.glVertex3d(p1.x,p1.y,p1.z);
				gl.glVertex3d(p2.x,p2.y,p2.z);
				gl.glVertex3d(p3.x,p3.y,p3.z);
			}
			gl.glEnd();
			gl.glLightModeli( GL2.GL_LIGHT_MODEL_TWO_SIDE, GL.GL_FALSE );
			gl.glEnable( GL.GL_CULL_FACE );
		}
	}

	private BooleanParameter useGravity = new BooleanParameter( "use gravity", true );
	private DoubleParameter gravity = new DoubleParameter( "gravity", 9.8, 0.01, 1000 );
	private DoubleParameter springStiffness = new DoubleParameter( "spring stiffness", 100, 0, 10000 );
	private DoubleParameter springDamping = new DoubleParameter( "spring damping", 0, 0, 50 );
	private DoubleParameter viscousDamping = new DoubleParameter( "viscous damping", 0, 0, 10 );
	private DoubleParameter restitution = new DoubleParameter( "r", 0, 0, 1 );
	private JTextArea comments = new JTextArea("enter comments in control panel");
	private DoubleParameter H = new DoubleParameter( "thickness", 0.5, 0.1, 10 );
	private IntParameter iterations = new IntParameter( "iterations", 10, 1, 100 );
	/** controls weather explicit or implicit integration is used */
	public BooleanParameter explicit = new BooleanParameter( "explicit", false );

	private IntParameter complexity = new IntParameter("cloth complexity", 15, 1, 50);
	private BooleanParameter pinTwoPoint = new BooleanParameter("two corner pinned", true);
	private BooleanParameter pinTop = new BooleanParameter("top edge pinned", false);
	private BooleanParameter drawGrid = new BooleanParameter( "draw spring grid", false );
	private BooleanParameter drawVertices = new BooleanParameter( "draw vertices", false );
	private BooleanParameter drawCloth = new BooleanParameter( "draw cloth", true );
	private BooleanParameter drawObj = new BooleanParameter( "draw objects - sphere", true);
	private DoubleParameter x1pos = new DoubleParameter( "x-coordinate of sphere 1", 1.8, -2, 2 );
	private DoubleParameter y1pos = new DoubleParameter( "y-coordinate of sphere 1", 0, -2, 2 );
	private DoubleParameter z1pos = new DoubleParameter( "z-coordinate of sphere 1", -.5, -2, 2 );
	private DoubleParameter radius1 = new DoubleParameter( "radius of sphere 1", .5, 0, 2 );
	private DoubleParameter x2pos = new DoubleParameter( "x-coordinate of sphere 2", 0, -2, 2 );
	private DoubleParameter y2pos = new DoubleParameter( "y-coordinate of sphere 2", 0, -2, 2 );
	private DoubleParameter z2pos = new DoubleParameter( "z-coordinate of sphere 2", 0, -2, 2 );
	private DoubleParameter radius2 = new DoubleParameter( "radius of sphere 2", 1, 0, 2 );
	private DoubleParameter radiusratio = new DoubleParameter( "radius display ratio", 0.9, .9, 1 );
	private BooleanParameter useSphere1 = new BooleanParameter( "use sphere 1" , true );
	private BooleanParameter useSphere2 = new BooleanParameter( "use sphere 2" , true );
	private BooleanParameter usePostStep = new BooleanParameter( "use post step contraints" , true );
	
	private BooleanParameter repulsion = new BooleanParameter( "apply repulsion impulses", true );
    
    private BooleanParameter collision = new BooleanParameter( "apply collision impulses", true );

    private IntParameter sphereSubdivision = new IntParameter( "sphere subdivision level", 2, 1, 5);

	/** parameters for controlling the position and size of an optional obstacle */

	private BooleanParameter flatShaded = new BooleanParameter( "flat shaded", false );

	@Override
	public JPanel getControls() {
		VerticalFlowPanel vfp = new VerticalFlowPanel();
		vfp.add( comments );
		vfp.add( useGravity.getControls() );
		vfp.add( gravity.getSliderControls(true) );
		vfp.add( springStiffness.getSliderControls(false) );
		vfp.add( springDamping.getSliderControls(false) );
		vfp.add( viscousDamping.getSliderControls(false) );
		vfp.add( restitution.getSliderControls(false) );
		vfp.add( H.getSliderControls(false));
		vfp.add( iterations.getSliderControls() );
		vfp.add( explicit.getControls() );
		vfp.add( pinTwoPoint.getControls());
		vfp.add( pinTop.getControls());
		vfp.add( complexity.getControls());
		vfp.add( drawGrid.getControls() );
		vfp.add( drawVertices.getControls() );
		vfp.add( drawCloth.getControls());
		vfp.add( drawObj.getControls());
		vfp.add( sphereSubdivision.getSliderControls());
		vfp.add( usePostStep.getControls());
		vfp.add( repulsion.getControls());
		vfp.add( collision.getControls());
		vfp.add( useSphere1.getControls());
		vfp.add( x1pos.getSliderControls(false) );
		vfp.add( y1pos.getSliderControls(false) );
		vfp.add( z1pos.getSliderControls(false) );
		vfp.add( radius1.getSliderControls(false) );
		vfp.add( flatShaded.getControls() );

		return vfp.getPanel();        
	}

	@Override
	public String toString() {
		String ret = "ZHIGUO ZHANG 260550226\n" +
				comments.getText() + "\n" +
				"particles = " + particles.size() + "\n";
		if ( explicit.getValue() ) {
			ret += "integrator = " + integrator.getName() + "\n";
		} else {
			ret += "integrator = Backward Euler\n";
		}
		ret += "k = " + springStiffness.getValue() + "\n" +
				"b = " + springDamping.getValue() + "\n" +
				"c = " + viscousDamping.getValue() +"\n" + 
				"time = " + time;
		return ret;
	}


	/**
	 * Fixes positions and velocities after a step to deal with collisions 
	 */
	public void postStep(int which) {
		Point3d pos = new Point3d();
		double radius = 0.0;
		if (which == 1) {
			pos.set(x1pos.getValue(), y1pos.getValue(), z1pos.getValue());
			radius = radius1.getValue();
		} else if (which == 2) {
			pos.set(x2pos.getValue(), y2pos.getValue(), z2pos.getValue());
			radius = radius2.getValue();
		} else return;

		for ( Particle p : particles ) {
			if ( p.pinned ) {
				p.v.set(0,0,0);
			}
		}
		// do sphere collisions
		double r = restitution.getValue();
//		for ( Triangle t : triangles) {
//			if (t.distance(pos.x, pos.y, pos.z) <= radius) {
//				Particle p = null;
//				for (int i = 0; i < 3; i++) {
//					if (i==0) p=t.p1;
//					else if (i==1) p=t.p2;
//					else p=t.p3;
//					Vector3d v = new Vector3d(p.p.x - pos.x, p.p.y - pos.y, p.p.z - pos.z);
//					v.normalize();
//					v.scale(radius);
//					p.p = new Point3d(pos.x + v.x, pos.y + v.y, pos.z + v.z);
//
//					if(p.v.length() > 0)
//					{
//						p.v.scale(r);
//						p.v.negate();
//					}
//					if(p.f.length() > 0)
//					{
//						p.f.set(0,0,0);
//					}
//				}
//				
//			}
//		}
		for ( Particle p : particles ) {            
			if(p.distance(pos.x, pos.y, pos.z) <= radius)
			{
				Vector3d v = new Vector3d(p.p.x - pos.x, p.p.y - pos.y, p.p.z - pos.z);
				v.normalize();
				v.scale(radius);
				p.p = new Point3d(pos.x + v.x, pos.y + v.y, pos.z + v.z);

				if(p.v.length() > 0)
				{
					p.v.scale(r);
					p.v.negate();
				}
				if(p.f.length() > 0)
				{
					p.f.set(0,0,0);
				}
			}
		}
	}

	/**
	 * Apply repulsion forces by particle vs. triangle test and boundary vs. boundary test
	 * @param h
	 */
	// TODO : Modify inverse mass and etc
	private void repulsion(double h) {
		// particle vs. triangle test
		for(Particle p: particles) {
			for(Triangle t: triangles) {
				if(p != t.p1 && p != t.p2 && p != t.p3) {
					Point3d p1 = t.p1.p;
					Point3d p2 = t.p2.p;
					Point3d p3 = t.p3.p;
					Vector3d v1 = new Vector3d();
					Vector3d v2 = new Vector3d();
					Vector3d normal = new Vector3d();

					v1.sub(p2,p1);
					v2.sub(p3,p1);
					normal.cross(v1,v2);
					normal.normalize();

					Vector3d v = new Vector3d();
					v.sub(p.p,p1);

					double dot = normal.dot(v);

					if(Math.abs(dot) < H.getValue()) {
						double l = v.dot(normal);
						normal.scale(l);

						Point3d project = new Point3d();
						project.sub(p.p, normal);

						Vector3d v0 = new Vector3d();
						v0.sub(project,p1);

						Vector3d normalA = new Vector3d();
						Vector3d normalB = new Vector3d();
						Vector3d normalC = new Vector3d();

						Vector3d v3 = new Vector3d();
						Vector3d v4 = new Vector3d();
						Vector3d v5 = new Vector3d();
						Vector3d v6 = new Vector3d();

						v3.sub(p3,p2);
						v4.sub(project,p2);
						v5.sub(p1,p3);
						v6.sub(project,p3);

						normalA.cross(v3, v4);
						normalB.cross(v5, v6);
						normalC.cross(v1, v0);

						normal.cross(v1,v2);

						double a = normal.dot(normalA)/normal.lengthSquared();
						double b = normal.dot(normalB)/normal.lengthSquared();
						double c = normal.dot(normalC)/normal.lengthSquared();

						normal.normalize();

						if(0 < a && a < 1 && 0 < b && b < 1 && 0 < c && c < 1)
						{
							Vector3d temp = new Vector3d();
							temp.sub(project,p.p);

							if(normal.dot(p.v) < 0) normal.negate();

							double d = H.getValue() - normal.dot(temp);

							Vector3d vRel = new Vector3d();
							Vector3d vTriangle = new Vector3d(a*t.p1.v.x + b*t.p2.v.x + c*t.p3.v.x, a*t.p1.v.y + b*t.p2.v.y + c*t.p3.v.y, a*t.p1.v.z + b*t.p2.v.z + c*t.p3.v.z);

							vRel.sub(p.v, vTriangle);

							if(vRel.dot(normal) >= 0.1*d/h)
							{
								double j = -Math.min(h*Spring.k*d, p.mass*(0.1*d/h - vRel.dot(normal)));
								double i = 2.0*j / (1.0+a*a+b*b+c*c);

								if (!p.pinned) p.v.add(new Vector3d(-i*normal.x/p.mass, -i*normal.y/p.mass, -i*normal.z/p.mass));
								if (!t.p1.pinned) t.p1.v.add(new Vector3d(i*a*normal.x/t.p1.mass, i*a*normal.y/t.p1.mass, i*a*normal.z/t.p1.mass));
								if (!t.p2.pinned) t.p2.v.add(new Vector3d(i*b*normal.x/t.p2.mass, i*b*normal.y/t.p2.mass, i*b*normal.z/t.p2.mass));
								if (!t.p3.pinned) t.p3.v.add(new Vector3d(i*c*normal.x/t.p3.mass, i*c*normal.y/t.p3.mass, i*c*normal.z/t.p3.mass));
							}
						}
					}
				}
			}
		}

		// edge vs edge test
		for (Spring s1 : springs) {
			for (Spring s2 : springs) {
				if ( s1 != s2 ) {
					// check if two springs intersect
					// http://paulbourke.net/geometry/pointlineplane/
					double EPS = 1e-6;
					Point3d p1, p2, p3, p4, pa, pb, p13, p43, p21;
					double mua, mub, d1343, d4321, d1321, d4343, d2121, numer, denom;
					p1 = s1.p1.p;
					p2 = s1.p2.p;
					p3 = s2.p1.p;
					p4 = s2.p2.p;
					p13 = new Point3d();
					p13.sub(p1, p3);
					p43 = new Point3d();
					p43.sub(p4, p3);
					if (Math.abs(p43.x) < EPS && Math.abs(p43.y) < EPS && Math.abs(p43.z) < EPS) continue;
					p21 = new Point3d();
					p21.sub(p2, p1);
					if (Math.abs(p21.x) < EPS && Math.abs(p21.y) < EPS && Math.abs(p21.z) < EPS) continue;
					d1343 = p13.x * p43.x + p13.y * p43.y + p13.z * p43.z;
					d4321 = p43.x * p21.x + p43.y * p21.y + p43.z * p21.z;
					d1321 = p13.x * p21.x + p13.y * p21.y + p13.z * p21.z;
					d4343 = p43.x * p43.x + p43.y * p43.y + p43.z * p43.z;
					d2121 = p21.x * p21.x + p21.y * p21.y + p21.z * p21.z;

					denom = d2121 * d4343 - d4321 * d4321;
					if (Math.abs(denom) < EPS) continue;
					numer = d1343 * d4321 - d1321 * d4343;

					mua = numer / denom;
					mub = (d1343 + d4321 * (mua)) / d4343;
					if (mua < 0 || mua > 1 || mub < 0 || mub > 1) continue;
					pa = new Point3d();
					pa.scaleAdd(mua, p21, p1);
					pb = new Point3d();
					pb.scaleAdd(mub, p43, p3);

					// TODO: finish collision
					Vector3d pab = new Vector3d();
					pab.sub(pb, pa);
					if (pab.length() < H.getValue()) {

						double d = H.getValue() - pab.length();

						Vector3d vRel = new Vector3d();
						Vector3d v1 = new Vector3d(mua * s1.p2.v.x + (1-mua) * s1.p1.v.x, mua * s1.p2.v.y + (1-mua) * s1.p1.v.y, mua * s1.p2.v.z + (1-mua) * s1.p1.v.z);
						Vector3d v2 = new Vector3d(mub * s2.p2.v.x + (1-mub) * s2.p1.v.x, mub * s2.p2.v.y + (1-mub) * s2.p1.v.y, mub * s2.p2.v.z + (1-mub) * s2.p1.v.z);

						vRel.sub(v1, v2);
						Vector3d normal = new Vector3d(pab);
						normal.normalize();
						if(vRel.dot(normal) >= 0.1*d/h)
						{
							double j = -Math.min(h*Spring.k*d, (s1.p1.mass + s1.p2.mass)*(0.1*d/h - vRel.dot(normal)));
							double i = 2.0*j / (mua*mua + (1-mua)*(1-mua) + mub*mub + (1-mub)*(1-mub));

							if (!s1.p1.pinned) s1.p1.v.add(new Vector3d(-(1-mua)*i*normal.x/s1.p1.mass,-(1-mua)*i*normal.y/s1.p1.mass,-(1-mua)*i*normal.z/s1.p1.mass));
							if (!s1.p2.pinned) s1.p2.v.add(new Vector3d(  -(mua)*i*normal.x/s1.p2.mass,  -(mua)*i*normal.y/s1.p2.mass,  -(mua)*i*normal.z/s1.p2.mass));
							if (!s2.p1.pinned) s2.p1.v.add(new Vector3d( (1-mub)*i*normal.x/s2.p1.mass, (1-mub)*i*normal.y/s2.p1.mass, (1-mub)*i*normal.z/s2.p1.mass));
							if (!s2.p2.pinned) s2.p2.v.add(new Vector3d(   (mub)*i*normal.x/s2.p2.mass,   (mub)*i*normal.y/s2.p2.mass,   (mub)*i*normal.z/s2.p2.mass));
						}
					}
				}
			}
		}
	}
	
	private List<Particle> allParticles = null;
    private List<Spring> allSprings = null;
    private List<Triangle> allTriangles = null;
    
	private boolean iterativeCheck(double h, int iter) {
		int MAX_ITERATION = iterations.getValue();
		double restitutiond = restitution.getValue();
		if (iter > MAX_ITERATION || iter == 0) return false; // give up
	    	boolean collision = false;
	    	
	    	allParticles = new LinkedList<Particle>(particles);
	    	allParticles.addAll(objParticles);
	    	
	    	allSprings = new LinkedList<Spring>(springs);
	    	allSprings.addAll(objSprings);
	    	
	    	allTriangles = new LinkedList<Triangle>(triangles);
	    	allTriangles.addAll(objTriangles);
	    	
	    	for(Particle p: allParticles) {
	    		for(Triangle t: allTriangles) {
	    			if(p != t.p1 && p != t.p2 && p != t.p3) {
					Point3d p1 = t.p1.p;
					Point3d p2 = t.p2.p;
					Point3d p3 = t.p3.p;
					Vector3d v1 = new Vector3d();
					Vector3d v2 = new Vector3d();
					Vector3d normal = new Vector3d();

					v1.sub(p2,p1);
					v2.sub(p3,p1);
					normal.cross(v1,v2);
					normal.normalize();

					Vector3d v = new Vector3d();
					v.sub(p.p,p1);

					double dot = normal.dot(v);
					if(Math.abs(dot) < 1e-6) {
						double l = v.dot(normal);
						normal.scale(l);

						Point3d project = new Point3d();
						project.sub(p.p, normal);

						Vector3d v0 = new Vector3d();
						v0.sub(project,p1);

						Vector3d normalA = new Vector3d();
						Vector3d normalB = new Vector3d();
						Vector3d normalC = new Vector3d();

						Vector3d v3 = new Vector3d();
						Vector3d v4 = new Vector3d();
						Vector3d v5 = new Vector3d();
						Vector3d v6 = new Vector3d();

						v3.sub(p3,p2);
						v4.sub(project,p2);
						v5.sub(p1,p3);
						v6.sub(project,p3);

						normalA.cross(v3, v4);
						normalB.cross(v5, v6);
						normalC.cross(v1, v0);

						normal.cross(v1,v2);

						double a = normal.dot(normalA)/normal.lengthSquared();
						double b = normal.dot(normalB)/normal.lengthSquared();
						double c = normal.dot(normalC)/normal.lengthSquared();

						normal.normalize();

						if(0 < a && a < 1 && 0 < b && b < 1 && 0 < c && c < 1) {
							Vector3d temp = new Vector3d();
							temp.sub(project,p.p);

							if(normal.dot(p.v) < 0) normal.negate();


							Vector3d vRel = new Vector3d();
							Vector3d vTriangle = new Vector3d(a*t.p1.v.x + b*t.p2.v.x + c*t.p3.v.x, a*t.p1.v.y + b*t.p2.v.y + c*t.p3.v.y, a*t.p1.v.z + b*t.p2.v.z + c*t.p3.v.z);

							vRel.sub(p.v, vTriangle);
							double vrel_neg = vRel.dot(normal);
							
							double i =  ( 1 + restitutiond ) * vrel_neg / (a * a / t.p1.mass + b * b / t.p2.mass+ c * c / t.p3.mass+ 1.0/p.mass);
									
							if (!p.pinned) p.v.add(new Vector3d(-i*normal.x/p.mass, -i*normal.y/p.mass, -i*normal.z/p.mass));
							if (!t.p1.pinned) t.p1.v.add(new Vector3d(i*a*normal.x/t.p1.mass, i*a*normal.y/t.p1.mass, i*a*normal.z/t.p1.mass));
							if (!t.p2.pinned) t.p2.v.add(new Vector3d(i*b*normal.x/t.p2.mass, i*b*normal.y/t.p2.mass, i*b*normal.z/t.p2.mass));
							if (!t.p3.pinned) t.p3.v.add(new Vector3d(i*c*normal.x/t.p3.mass, i*c*normal.y/t.p3.mass, i*c*normal.z/t.p3.mass));
							collision = true;
						}
					}
	    			}
			}
	    	}
	    	
	    	// edge vs edge test
		for (Spring s1 : allSprings) {
			for (Spring s2 : allSprings) {
				if ( s1 != s2 ) {
					// check if two springs intersect
					// http://paulbourke.net/geometry/pointlineplane/
					double EPS = 1e-6;
					Point3d p1, p2, p3, p4, pa, pb, p13, p43, p21;
					double mua, mub, d1343, d4321, d1321, d4343, d2121, numer, denom;
					p1 = s1.p1.p;
					p2 = s1.p2.p;
					p3 = s2.p1.p;
					p4 = s2.p2.p;
					p13 = new Point3d();
					p13.sub(p1, p3);
					p43 = new Point3d();
					p43.sub(p4, p3);
					if (Math.abs(p43.x) < EPS && Math.abs(p43.y) < EPS && Math.abs(p43.z) < EPS) continue;
					p21 = new Point3d();
					p21.sub(p2, p1);
					if (Math.abs(p21.x) < EPS && Math.abs(p21.y) < EPS && Math.abs(p21.z) < EPS) continue;
					d1343 = p13.x * p43.x + p13.y * p43.y + p13.z * p43.z;
					d4321 = p43.x * p21.x + p43.y * p21.y + p43.z * p21.z;
					d1321 = p13.x * p21.x + p13.y * p21.y + p13.z * p21.z;
					d4343 = p43.x * p43.x + p43.y * p43.y + p43.z * p43.z;
					d2121 = p21.x * p21.x + p21.y * p21.y + p21.z * p21.z;

					denom = d2121 * d4343 - d4321 * d4321;
					if (Math.abs(denom) < EPS) continue;
					numer = d1343 * d4321 - d1321 * d4343;

					mua = numer / denom;
					mub = (d1343 + d4321 * (mua)) / d4343;
					if (mua < 0 || mua > 1 || mub < 0 || mub > 1) continue;
					pa = new Point3d();
					pa.scaleAdd(mua, p21, p1);
					pb = new Point3d();
					pb.scaleAdd(mub, p43, p3);

					// TODO: finish collision
					Vector3d pab = new Vector3d();
					pab.sub(pb, pa);
					if (pab.length() < EPS) {
						Vector3d vRel = new Vector3d();
						Vector3d v1 = new Vector3d(mua * s1.p2.v.x + (1-mua) * s1.p1.v.x, mua * s1.p2.v.y + (1-mua) * s1.p1.v.y, mua * s1.p2.v.z + (1-mua) * s1.p1.v.z);
						Vector3d v2 = new Vector3d(mub * s2.p2.v.x + (1-mub) * s2.p1.v.x, mub * s2.p2.v.y + (1-mub) * s2.p1.v.y, mub * s2.p2.v.z + (1-mub) * s2.p1.v.z);

						vRel.sub(v1, v2);
						Vector3d normal = new Vector3d(pab);
						normal.normalize();
						double vrel_neg = vRel.dot(normal);
						double i =  ( 1 + restitutiond ) * vrel_neg / (mua*mua/s1.p2.mass + (1-mua)*(1-mua)/s1.p1.mass + mub*mub/s2.p2.mass + (1-mub)*(1-mub)/s2.p1.mass);

						if (!s1.p1.pinned) s1.p1.v.add(new Vector3d(-(1-mua)*i*normal.x/s1.p1.mass,-(1-mua)*i*normal.y/s1.p1.mass,-(1-mua)*i*normal.z/s1.p1.mass));
						if (!s1.p2.pinned) s1.p2.v.add(new Vector3d(  -(mua)*i*normal.x/s1.p2.mass,  -(mua)*i*normal.y/s1.p2.mass,  -(mua)*i*normal.z/s1.p2.mass));
						if (!s2.p1.pinned) s2.p1.v.add(new Vector3d( (1-mub)*i*normal.x/s2.p1.mass, (1-mub)*i*normal.y/s2.p1.mass, (1-mub)*i*normal.z/s2.p1.mass));
						if (!s2.p2.pinned) s2.p2.v.add(new Vector3d(   (mub)*i*normal.x/s2.p2.mass,   (mub)*i*normal.y/s2.p2.mass,   (mub)*i*normal.z/s2.p2.mass));
						collision = true;
					}
				}
			}
		}
		
	    	if (collision) return iterativeCheck(h, iter+1); // Might be more collision!
	    	else return true; // No more collision found!
    }
}
