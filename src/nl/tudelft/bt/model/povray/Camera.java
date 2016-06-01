package nl.tudelft.bt.model.povray;

import java.io.Serializable;

/**
 * Camera object
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class Camera  implements Serializable {
	private VectorProperty location;
	private VectorProperty up;
	private VectorProperty right;
	private VectorProperty look_at;
	private float angle;

	public Camera() {
		location = new VectorProperty("location");
		up = new VectorProperty("up");
		right = new VectorProperty("right");
		look_at = new VectorProperty("look_at");
	}

	/**
	 * @param fs
	 */
	public void setLocation(float x, float y, float z) {
		location.setValues(x, y, z);
	}

	public void setUp(float x, float y, float z) {
		up.setValues(x, y, z);
	}

	public void setRight(float x, float y, float z) {
		right.setValues(x, y, z);
	}

	/**
	 * @param fs
	 */
	public void setLook_at(float x, float y, float z) {
		look_at.setValues(x, y, z);
	}

	/**
	 * @param f
	 */
	public void setAngle(float f) {
		angle = f;
	}

	public String toString() {
		return "camera {\n"
			+ "\t"
			+ location
			+ "\n"
			+ "\t "
			+ up
			+ "\n"
			+ "\t "
			+ right
			+ "\n"
			+ "\t "
			+ look_at
			+ "\n"
			+ "\tangle "
			+ angle
			+ "\n"
			+ "}\n";
	}

}
