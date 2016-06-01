/*
 * Created on Jun 18, 2003
 *
 * To change the template for this generated file go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */
package nl.tudelft.bt.model.povray;
import java.io.*;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.multigrid.boundary_conditions.BoundaryConditions;
/**
 * Creates a full 3D scene for povRay from the bacteria container
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class Povray3DScene implements Serializable {
	private static Model _model = Model.model();
	final private static String INCHEADER = "sceneheader.inc";
	final private static String INCFOOTER = "scenefooter.inc";
	private Camera _camera;
	private Background _background;
	private LightSource[] _lightSource;
	private Biofilm3D _biofilm;
	private float _x;
	private float _y;
	private float _z;
	private static float _scalling;
	public Povray3DScene() {
		_scalling = _model.systemSize.y;
		_x = _model.systemSize.x / _scalling;
		_y = _model.systemSize.y / _scalling;
		_z = _model.systemSize.z / _scalling;
		initializeCameraBackgroundAndLight();
	}
	/**
	 * Constructs a povray scene for a system size of x by y by z and n cells
	 * 
	 * @param x
	 * @param y
	 * @param z
	 * @param n
	 */
	public Povray3DScene(float x, float y, float z, int n) {
		_x = x;
		_y = y;
		_z = z;
		initializeCameraBackgroundAndLight();
		_biofilm = new Biofilm3D(this, n);
	}
	/**
	 * Make top perspective
	 */
	public void setTopPerspective() {
		_camera.setLocation(0, _y * 1.56f, 0);
		_lightSource[1].setLocation(_z, 0, 0);
	}
	/**
	 * Make side perspective
	 */
	public void setSidePerspective() {
		_camera.setLocation(0, _y, _x * 1.17f);
		_lightSource[1].setLocation(_z, 0, 0);
		_camera.setLook_at(0, -_y * 0.5f, 0);
	}
	/**
	 * Make angle perspective
	 */
	public void setAnglePerspective() {
		_camera.setLocation(_z, _y, _x * 2.5f);
		_lightSource[1].setLocation(_z, 0, 0);
		_camera.setLook_at(0, 0, 0);
	}
	private void initializeCameraBackgroundAndLight() {
		_camera = new Camera();
		if (Model.model().getDimensionality() == 3)
			_camera.setLocation(0, _y * 2f, _x * 2f);
		else
			_camera.setLocation(0, 0, _x * 1.7f);
		_camera.setUp(0, 1, 0);
		_camera.setRight(-1.33f, 0, 0);
		_camera.setLook_at(0, 0, 0);
		_camera.setAngle(40);
		_background = new Background();
		_background.setColor(1f, 1f, 1f);
		_lightSource = new LightSource[2];
		_lightSource[0] = new LightSource();
		_lightSource[0].setLocation(_z, _y, _x);
		_lightSource[0].setColor(1f, 1f, 1f);
		_lightSource[1] = new LightSource();
		_lightSource[1].setLocation(-_z, _y, _x);
		_lightSource[1].setColor(1f, 1f, 1f);
	}
	/**
	 * Add a cell to the biofilm, to be used with scenes constructed using the
	 * Povray3DScene(float x, float y, float z, int n) constructor
	 * 
	 * @param x
	 *            coordinate
	 * @param y
	 *            coordinate
	 * @param z
	 *            coordinate
	 * @param rad
	 *            radius
	 * @param r
	 *            red
	 * @param g
	 *            green
	 * @param b
	 *            blue
	 */
	public void addCell(float x, float y, float z, float rad, int r, int g,
			int b) {
		_biofilm.addCell(x, y, z, rad, r, g, b);
	}
	/**
	 * Write the complete model state to a file (including camera, lights and
	 * boxes)
	 * 
	 * @param fn
	 * @throws IOException
	 */
	public void writeModelStateToPovrayFileFull(String fn) throws IOException {
		// update the biofilm particles
		_biofilm = new Biofilm3D(this);
		java.io.File sysData = new java.io.File(fn);
		FileWriter fr = new FileWriter(sysData);
		fr.write(_camera.toString());
		fr.write(_background.toString());
		fr.write(_lightSource[0].toString());
		fr.write(_lightSource[1].toString());
		_biofilm.modelStateToFile(fr);
		fr.close();
	}
	/**
	 * Write the include files. For use with writeModelStateToPovrayFile
	 * 
	 * @param dir
	 *            directory to write the include file
	 * @throws IOException
	 */
	public File[] writePovrayIncFiles(String dir) throws IOException {
		_biofilm = new Biofilm3D(this);
		//header include file
		java.io.File header = new java.io.File(dir + INCHEADER);
		FileWriter fr = new FileWriter(header);
		fr.write(_camera.toString());
		fr.write(_background.toString());
		fr.write(_lightSource[0].toString());
		fr.write(_lightSource[1].toString());
		_biofilm.biofilmHeaderToFile(fr);
		fr.close();
		//footer include file
		java.io.File footer = new java.io.File(dir + INCFOOTER);
		fr = new FileWriter(footer);
		_biofilm.biofilmFooterToFile(fr);
		fr.close();
		//prepare file array to return
		File[] incs = {header, footer};
		return incs;
	}
	/**
	 * Write the present state using include files for camera, background lights
	 * and solid surface. Using this method instead of
	 * writeModelStateToPovrayFileFull, which includes the full scene
	 * information in each .pov file created at an iteration, allows changing
	 * scene properties after the simulation
	 * 
	 * @param fn
	 * @return the file written
	 * @throws IOException
	 */
	public File writeModelStateToPovrayFile(String fn) throws IOException {
		// update the biofilm particles
		_biofilm = new Biofilm3D(this);
		File sysData = new File(fn);
		FileWriter fr = new FileWriter(sysData);
		fr.write("#include \"" + INCHEADER + "\"\n");
		_biofilm.particlesToFile(fr);
		fr.write("#include \"" + INCFOOTER + "\"\n");
		fr.close();
		return sysData;
	}
	/**
	 * @return
	 */
	protected float getX() {
		return _x;
	}
	/**
	 * @return
	 */
	protected float getY() {
		return _y;
	}
	/**
	 * @return
	 */
	protected float getZ() {
		return _z;
	}
	/**
	 * @return Returns the _scalling.
	 */
	public static float getScalling() {
		return _scalling;
	}
	
	/**
	 * @param os
	 * @throws IOException
	 */
	public static void serializeStaticState(ObjectOutputStream os)
			throws IOException {
		os.writeFloat(_scalling);
	}

	/**
	 * @param os
	 * @throws IOException
	 * @throws ClassNotFoundException
	 */
	public static void deserializeStaticState(ObjectInputStream os)
			throws IOException, ClassNotFoundException {
		_scalling = os.readFloat();
	}

}