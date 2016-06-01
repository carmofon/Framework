package nl.tudelft.bt.model.work.drift;

import java.awt.Color;
import java.io.Serializable;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.BiomassSpecies.Composition;
import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.particlebased.BiomassParticleContainer;
import nl.tudelft.bt.model.util.ColorMaps;

public class BiomassWithNeutralColors extends BiomassSpecies {

	public BiomassWithNeutralColors(String name, ParticulateSpecies[] species,
			float[] fractionalCompositionInVolume)
			throws NonMatchingNumberException {
		super(name, species, fractionalCompositionInVolume);
		// TODO Auto-generated constructor stub
	}

	/**
	 * May be overriden to call the instance specific inner class Composition to
	 * allow more flexibility in extending BiomassSpecies
	 * 
	 * @param s
	 * @param masses
	 * @return the instance of composition
	 */
	protected Composition createComposition(BiomassSpecies s, float[] masses) {
		return new Composition2(this, masses);
	}

	/**
	 * Returns a new particle with random volume of 80 to 90 %of the value
	 * defined by the maximum radius
	 * 
	 * @return a volumetric composition composition
	 * @throws NonMatchingNumberException
	 *             never occurs, as number consistency is checked in the
	 *             constructor
	 */
	public Composition newParticleComposition()
			throws NonMatchingNumberException {
		Composition c = super.newParticleComposition();
		((Composition2)c).setRandomColor();
		return c;
	}

	/**
	 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
	 */
	public class Composition2 extends Composition {
		private Color _color;

		public Composition2(BiomassSpecies s, float[] masses)
				throws NonMatchingNumberException {
			super(s, masses);
			// TODO Auto-generated constructor stub
		}

		public Composition2(Composition c) {
			super(c);
			// TODO Auto-generated constructor stub
		}

		@Override
		public Color getColorCore() {
			// TODO Auto-generated method stub
			return _color;
		}

		@Override
		protected Object clone() throws CloneNotSupportedException {
			// TODO Auto-generated method stub
			return super.clone();
		}
		
		private void setRandomColor() {
			_color = ColorMaps.getJetColor(Model.model().getRandom());
		}
	}

}
