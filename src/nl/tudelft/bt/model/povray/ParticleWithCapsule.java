/*
 * Created on Jun 18, 2003
 * 
 * To change the template for this generated file go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */
package nl.tudelft.bt.model.povray;

import java.awt.Color;
import java.io.Serializable;

import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.util.ColorMaps;

/**
 * Sphere object for povray scene
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ParticleWithCapsule  implements Serializable {
	private VectorProperty center;

	private float _coreRadius;

	private float _capsuleRadius;

	private VectorProperty _colorCore;

	private VectorProperty _colorCapsule;

	private boolean _hasCapsule;

	public ParticleWithCapsule() {
		center = new VectorProperty("");
		_colorCore = new VectorProperty("color rgb");
	}

	public ParticleWithCapsule(BiomassParticle p) {
		center = new VectorProperty("");
		setCenter(p.getCenter());
		_colorCore = new VectorProperty("color rgb");
		setColorCore(p.getColorCore());
		setCoreRadius(p.getCoreRadius());
		_hasCapsule = p.hasCapsule();
		if (_hasCapsule) {
			_capsuleRadius = p.getRadius()/Povray3DScene.getScalling();
			_colorCapsule = new VectorProperty("rgbf");
			setColorCapsule(p.getColorCapsule());
		}
	}

	/**
	 * @param color
	 */
	public void setColorCore(Color c) {
		_colorCore.setValues(((float) c.getRed()) / 255,
				((float) c.getGreen()) / 255, ((float) c.getBlue()) / 255);
	}

	/**
	 * For now sets capsule to gray
	 * 
	 * @param fs
	 */
	public void setColorCapsule(Color c) {
		float r = ColorMaps.brightenValue(((float) c.getRed()) / 255, 0.5f);
		float g = ColorMaps.brightenValue(((float) c.getGreen()) / 255, 0.5f);
		float b = ColorMaps.brightenValue(((float) c.getBlue()) / 255, 0.5f);
		_colorCapsule.setValues(r, g, b, 0.999f);
	}

	/**
	 * @param fs
	 */
	public void setCenter(ContinuousCoordinate c) {
		float s = Povray3DScene.getScalling();
		center.setValues(c.x/s, c.y/s, c.z/s);
	}

	/**
	 * @param fs
	 */
	public void setCoreRadius(float fs) {
		_coreRadius = fs/Povray3DScene.getScalling();
	}

	public String toString() {
		String core = "sphere {\n" + "\t " + center + "\n" + "\t "
				+ _coreRadius + "\n" + "\t pigment { " + _colorCore + " }\n"
				+ "}\n";
		if (_hasCapsule) {
			String capsule = "sphere {\n" + "\t " + center + "\n" + "\t "
					+ _capsuleRadius + "\n" + "\t pigment { " + _colorCapsule
					+ " }\n" + "}\n";
			return core + capsule;
		}
		return core;
	}
}