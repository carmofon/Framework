package nl.tudelft.bt.model.work.drift;

import java.awt.Color;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

public class BiomassParticleInheritRandomColor extends BiomassParticle {

	public BiomassParticleInheritRandomColor(BiomassSpecies s) {
		super(s);
	}

	/**
	 * Returns the color for this bacterium, based on its growth rate
	 * 
	 * @return the color species of the bacteria (or _color if overidden)
	 */
	public Color getColorCore() {
		return Color.magenta;
	}

}
