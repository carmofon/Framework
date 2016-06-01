/*
 * Created on 8-jan-2004 by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model;

import java.awt.Color;
import java.io.Serializable;
import java.util.*;

import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.MultigridVariable;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.particlebased.BiomassParticleContainer;
import nl.tudelft.bt.model.reaction.Reaction;
import nl.tudelft.bt.model.util.ColorMaps;
import nl.tudelft.bt.model.util.ExtraMath;

/**
 * Implements the biomass species. Each biomass species will be composed by a
 * set of fixed species. For example, a biomass species may contain active
 * heterotroph biomass, inert heterotroph biomass EPS, etc.
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 * @see ParticulateSpecies
 */
public class BiomassSpecies implements Serializable {
	private String _name;

	protected ParticulateSpecies[] _particulates;

	protected float[] _newParticleFractionalVolume;

	// The following references are implemented for definig capsule surrounding
	// the biomass core
	private ParticulateSpecies _activeMass;

	private ParticulateSpecies _inertMass;

	private ParticulateSpecies[] _eps;

	private Reaction[] _reactions;

	protected boolean _getColorFromGrowth;

	private int _createdParticles;

	private boolean _overideColor = false;

	private SoluteSpecies _inducer;

	private float _induceThreshold;

	private Color _inducedColor;

	private Color _uninducedColor;

	private int _shovingHierarchy = 1;

	private static final float GAMMA = 0.4f;

	// precision of integration: t is subdivided into NINT subintervals
	public static final int NINT = 50;

	/**
	 * Create new biomass species using an array of fixed species
	 * 
	 * @param name
	 *            name of biomass species
	 * @param species
	 *            array of composing fixed species
	 * @param fractionalCompositionInVolume
	 *            array of fractional composition, sum of all entries must be 1
	 *            and number of entries must be the same as lenght of fixed
	 *            species array
	 * @throws NonMatchingNumberException
	 */
	public BiomassSpecies(String name, ParticulateSpecies[] species,
			float[] fractionalCompositionInVolume)
			throws NonMatchingNumberException {
		// check that no repetition exists
		assertThatNoRepetitionExists(species);
		// check length are consistent
		if (fractionalCompositionInVolume.length != species.length)
			throw new NonMatchingNumberException("Length of masses array"
					+ " inconsistent with number of fixed species");
		// check that sum of compositions is less than 1
		float total = 0;
		for (int i = 0; i < fractionalCompositionInVolume.length; i++) {
			total += fractionalCompositionInVolume[i];
		}
		if (total > 1)
			throw new NonMatchingNumberException("Sum of composition is"
					+ " larger than 1");
		// initialize the attributes
		_name = name;
		_particulates = species;
		_newParticleFractionalVolume = fractionalCompositionInVolume;
		_createdParticles = 0;
	}

	public void setInducibleColor(SoluteSpecies inducer, float threshold,
			Color inducedColor, Color unInducedColor) {
		_overideColor = true;
		_inducer = inducer;
		_induceThreshold = threshold;
		_inducedColor = inducedColor;
		_uninducedColor = unInducedColor;
	}

	/**
	 * @return Returns the _name.
	 */
	public String getName() {
		return _name;
	}

	/**
	 * @return a new biomass particle of this species
	 */
	public BiomassParticle createBiomassParticle() {
		return new BiomassParticle(this);
	}

	/**
	 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
	 */
	public class Composition implements Cloneable, Serializable {
		private BiomassSpecies _biomassSp;

		protected ParticulateSpecies[] _particulateSp;

		protected float[] _masses;

		private float[] _growthRates;

		private float _volume;

		private float _volumetricGrowthRate;

		private float _growthToMaximum;

		/**
		 * Creates an instance of composition according to a certain biomass
		 * species
		 * 
		 * @param s
		 * @param masses
		 * @throws NonMatchingNumberException
		 */
		public Composition(BiomassSpecies s, float[] masses)
				throws NonMatchingNumberException {
			int n = s._particulates.length;
			if (masses.length != n)
				throw new NonMatchingNumberException("Length of masses array ("
						+ masses.length
						+ ") inconsistent with number of fixed species (" + n
						+ ")");
			_biomassSp = s;
			_particulateSp = s._particulates;
			_masses = masses;
			_growthRates = new float[n];
		}

		/**
		 * Create a new instance as a copy of composition were all masses are
		 * set to 0
		 * 
		 * @param c
		 *            composition to copy
		 */
		public Composition(Composition c) {
			int n = c._biomassSp._particulates.length;
			_biomassSp = c._biomassSp;
			_particulateSp = _biomassSp._particulates;
			_masses = new float[n];
			_growthRates = new float[n];
			_growthToMaximum = c._growthToMaximum;
		}

		/**
		 * Compute the radius of a spherical biomass particle with this
		 * composition (including core and capsule) and also update the private
		 * attribute _volume;
		 * 
		 * @return radius of spherical particle
		 */
		public float computeRadius() {
			_volume = 0.0f;
			for (int i = 0; i < _particulateSp.length; i++) {
				_volume += _masses[i] / _particulateSp[i].getDensity();
			}
			return computeRadius(_volume);
		}

		/**
		 * Determine the radius of sphericla particle core (all masses except
		 * eps)
		 * 
		 * @return radius of spherical particle core
		 */
		public float getCoreRadius() {
			float volume = 0.0f;
			for (int i = 0; i < _particulateSp.length; i++) {
				// add volume corresponding to _particulateSp[i] only if it is
				// not contained in _eps array
				if (!containedInArray(_particulateSp[i], _eps))
					volume += _masses[i] / _particulateSp[i].getDensity();
			}
			return computeRadius(volume);
		}

		/**
		 * Compute the radius in 2D or 3D cases
		 * 
		 * @param v
		 * @return
		 */
		private float computeRadius(float v) {
			if (Model.model().getDimensionality() == 3)
				// 3D case
				// radius is computed from a sphere
				return ExtraMath.radiusOfASphere(v);
			else
				// 2D case
				// radius is computed from a cylinder
				return ExtraMath.radiusOfACilinder(v, Model.model()
						.get2DSystem3rdDimension());
		}

		/**
		 * Compute the growth rates for each of the fixed species and update the
		 * precomputed growth rates array
		 * 
		 * @param c
		 *            position of center of biomass particle
		 */
		public void computeGrowthRates(ContinuousCoordinate c) {
			// check if the reactions array is initialized. If not,
			// initialize it
			if (_biomassSp._reactions == null) {
				_biomassSp.initializeReactionsArray();
			}
			// precompute the rates for all reactions
			for (int i = 0; i < _biomassSp._reactions.length; i++) {
				MultigridVariable.setCurrentLocation(c);
				_biomassSp._reactions[i]
						.computeMassGrowthRateAndAddToGlobal(this);
			}
			// iterate through the reactions involved with this biomass
			// species and call getMassRate which performes all the
			// stoichiometries
			for (int i = 0; i < _particulateSp.length; i++) {
				_growthRates[i] = _particulateSp[i].getMassRate();
			}
		}

		/**
		 * Get the value for the volumetric growth rate
		 * 
		 * @return volumetric growth rate [um3/h]
		 */
		public float getVolumetricGrowthRate() {
			return _volumetricGrowthRate;
		}

		/**
		 * Computes the time constraint imposed by a decreasing rate. The time
		 * constraint imposes that time step should not be higher than the value
		 * to take the concentration to 0
		 * 
		 * @return the time constraint for particle with this composition
		 */
		public float getMaximumTimeConstraint() {
			float t = Float.POSITIVE_INFINITY;
			_volumetricGrowthRate = 0.0f;
			float volume = 0.0f;
			// determine the full volume and the total growth rate
			for (int i = 0; i < _particulateSp.length; i++) {
				// volume
				volume += _masses[i] / _particulateSp[i].getDensity();
				// growth rate sum
				_volumetricGrowthRate += _growthRates[i]
						/ _particulateSp[i].getDensity();
				// TEST: remove the constraint for particle shrinking
				// //BEGINIG
				// // compute the consumption time
				// if (_growthRates[i] < 0) {
				// //consumption in one time step must not be greater than
				// // the total mass of the species in one particle
				// t = Math.min(t, - 2f * _masses[i] / _growthRates[i]);
				// }
				// //END
			}
			// if volume change is positive, constraint is also that volume
			// should not double at each iteration
			// 0.69 (= ln 2) ensures this condition also for exponential growth
			// 0.5 (< ln 2) works better
			if (_volumetricGrowthRate > 0)
				return Math.min(t, volume / _volumetricGrowthRate * 0.5f);
			// otherwise, only the non-null mass restraint is considered
			return t;
		}

		/**
		 * Update the composition by integrating the growth rates along a time
		 * interval t
		 * 
		 * @param t
		 *            time step to integrate growth
		 * @param c
		 *            position of the center
		 * @return the updated value for the mass
		 */
		public float grow(float t, ContinuousCoordinate c) {
			// final int NINT = 1;
			float m = 0;
			_growthToMaximum = 1;
			// Temporary variable must be used to ensure consistence in
			// computation of all masses
			float[] newMasses = new float[_masses.length];
			// size of sub-interval
			float tsubint = t / NINT;
			for (int j = 0; j < NINT; j++) {
				for (int i = 0; i < _particulateSp.length; i++) {
					// compute the new mass
					// growth rate must have been previously computed
					newMasses[i] = _masses[i] + tsubint * _growthRates[i];
					// if lower than 0, truncate to 0
					if (newMasses[i] < 0)
						newMasses[i] = 0;
					// add mass change to total produced mass
					m += newMasses[i] - _masses[i];
					_masses[i] = newMasses[i];
				}
				// update growth rates with new particle biomass
				computeGrowthRates(c);
			}
			// update _growthToMaximum
			for (int i = 0; i < _particulateSp.length; i++) {
				_growthToMaximum *= _particulateSp[i].getSpecificRate(c)
						/ _particulateSp[i]
								.getCurrentMaximumSpecificGrowthRate();
			}
			return m;
		}

		/**
		 * Divide the masses asymmetrically and return a new instance of
		 * Composition with the remainder masses
		 * 
		 * @param f
		 *            fraction for division
		 * @return a new instance of composition containing the remainder of
		 *         masses
		 */
		public Composition divideMasses(float f) {
			Composition c = new Composition(this);
			for (int i = 0; i < _masses.length; i++) {
				c._masses[i] = _masses[i] * (1 - f);
				_masses[i] *= f;
			}
			return c;
		}

		/**
		 * Remove the eps capsule, creating a new capsule with EPS only
		 * 
		 * @return the composition of new particle
		 */
		public Composition removeCapsule() {
			Composition c = new Composition(this);
			// set mass of all eps components in new particle to match those
			// of old particle
			for (int i = 0; i < _biomassSp._eps.length; i++) {
				ParticulateSpecies epsComponent = _biomassSp._eps[i];
				// set mass in eps only particle
				c.setMass(epsComponent, getSpeciesMass(epsComponent));
				// set to 0 on this particle
				this.setMass(epsComponent, 0);
			}
			return c;
		}

		/**
		 * Sets the value of mass of a given fixed species
		 * 
		 * @param s
		 *            species to change value
		 * @param m
		 *            new mass value
		 * @throws ModelException
		 *             if species is not contained in the composition
		 */
		public void setMass(ParticulateSpecies s, float m)
				throws ModelRuntimeException {
			for (int i = 0; i < _particulateSp.length; i++) {
				if (_particulateSp[i] == s) {
					_masses[i] = m;
					return;
				}
			}
			throw new ModelRuntimeException("Trying to edit value of a fixed"
					+ " species not contained in this biomass species");
		}

		/**
		 * Add contribution of this composition to discrete data matrices
		 * 
		 * @param c
		 *            center of particle
		 */
		public void addContributionToDiscreteData(ContinuousCoordinate c) {
			for (int i = 0; i < _particulateSp.length; i++) {
				_particulateSp[i].addContributionToDiscreteData(c, _masses[i]);
			}
		}

		/**
		 * Return a color for the core of a particle which is the weighed
		 * average of the colors of the composing particulates, except eps,
		 * which composes the capsule, and inert, which causes the color to turn
		 * gray
		 * 
		 * @return color for particle with this composition
		 */
		public Color getColorCore(ContinuousCoordinate center) {
			boolean STEP = false;
			// check for color override at the biomass species level
			if (_overideColor) {
				float s = _inducer.getValueAt(center);
				// TODO switch this to have step like working of color
				if (STEP)
					return (s > _induceThreshold ? _inducedColor
							: _uninducedColor);
				// produces intermediate between colors
				float f = s / (s + _induceThreshold);
				float iRed = ((float) (_inducedColor.getRed())) / 255f;
				float iGreen = ((float) (_inducedColor.getGreen())) / 255f;
				float iBlue = ((float) (_inducedColor.getBlue())) / 255f;
				float uRed = ((float) (_uninducedColor.getRed())) / 255f;
				float uGreen = ((float) (_uninducedColor.getGreen())) / 255f;
				float uBlue = ((float) (_uninducedColor.getBlue())) / 255f;
				return new Color(iRed * f + uRed * (1 - f), iGreen * f + uGreen
						* (1 - f), iBlue * f + uBlue * (1 - f));

			}
			//
			if (_biomassSp._getColorFromGrowth)
				return ColorMaps.getFullJetColor(getRelativeGrowth(), 0.1f);
			//return ColorMaps.getRedscaleColor(getRelativeGrowth());
			float h = 0;
			float s = 0;
			float b = 0;
			// compute total mass of core of particle, i.e. excluding the eps
			// capsule and inert mass
			float totalMassInCore = 0;
			for (int i = 0; i < _masses.length; i++) {
				// check if species is not part of EPS nor inert
				if (!containedInArray(_particulateSp[i], _biomassSp._eps)
						& (_particulateSp[i] != _biomassSp._inertMass)) {
					totalMassInCore += _masses[i];
				}
			}
			// interpolate color from weighed average of composing particulate
			// species
			for (int i = 0; i < _masses.length; i++) {
				// check if species is not part of EPS nor inert
				if (!containedInArray(_particulateSp[i], _biomassSp._eps)
						& (_particulateSp[i] != _biomassSp._inertMass)) {
					float f = _masses[i] / totalMassInCore;
					Color c = _particulateSp[i].getColor();
					float[] hsb = new float[3];
					Color.RGBtoHSB(c.getRed(), c.getGreen(), c.getBlue(), hsb);
					h += hsb[0] * f;
					s += hsb[1] * f;
					b += hsb[2] * f;
				}
			}
			// Return right now if color is black
			if ((h == 0) & (s == 0) & (b == 0))
				return Color.black;
			// inert fraction causes the color to turn gray, without changing
			// the hue:
			float inertMass = getSpeciesMass(_biomassSp._inertMass);
			// decrase saturation using a gamma correction to emphasize color
			try {
				s = ExtraMath.gammaCorrection(1 - inertMass
						/ (totalMassInCore + inertMass), GAMMA);
			} catch (ModelException e) {
				e.printStackTrace();
				throw new ModelRuntimeException(e.toString());
			}
			// decrease brightness to a minimum of 75 %
			b = 1 - inertMass / (totalMassInCore + inertMass) * .15f;
			return Color.getHSBColor(h, s, b);
		}

		/**
		 * Retrun the color of the capsule (for instance EPS) around the
		 * particle.
		 * 
		 * @return color for capsule
		 */
		public Color getColorCapsule() {
			float h = 0;
			float s = 0;
			float b = 0;
			// compute total mass of capsule of particle
			float totalMassInCapsule = 0;
			for (int i = 0; i < _masses.length; i++) {
				// check that species is part of EPS
				if (containedInArray(_particulateSp[i], _biomassSp._eps)) {
					totalMassInCapsule += _masses[i];
				}
			}
			// interpolate color from weighed average of composing particulate
			// species
			for (int i = 0; i < _masses.length; i++) {
				// check if species is not part of EPS nor inert
				if (containedInArray(_particulateSp[i], _biomassSp._eps)) {
					float f = _masses[i] / totalMassInCapsule;
					Color c = _particulateSp[i].getColor();
					float[] hsb = new float[3];
					Color.RGBtoHSB(c.getRed(), c.getGreen(), c.getBlue(), hsb);
					h += hsb[0] * f;
					s += hsb[1] * f;
					b += hsb[2] * f;
				}
			}
			return Color.getHSBColor(h, s, b);
		}

		/**
		 * @param s
		 *            fixed species to get mass
		 * @return the mass on species s, 0 if s is not in the composition
		 */
		public float getSpeciesMass(ParticulateSpecies s) {
			for (int i = 0; i < _particulateSp.length; i++) {
				if (s == _particulateSp[i])
					return _masses[i];
			}
			return 0;
		}

		/**
		 * @param s
		 * @return true is species s is part of the biomass composition
		 */
		public boolean speciesPartOfComposition(ParticulateSpecies s) {
			for (int i = 0; i < _particulateSp.length; i++) {
				if (s == _particulateSp[i])
					return true;
			}
			return false;
		}

		/**
		 * @return the total mass of a particle
		 */
		public float getTotalMass() {
			float m = 0;
			for (int i = 0; i < _masses.length; i++) {
				m += _masses[i];
			}
			return m;
		}

		/**
		 * @return the total mass of EPS components in this particle
		 */
		public float getEpsMass() {
			float m = 0;
			for (int i = 0; i < _biomassSp._eps.length; i++) {
				m += getSpeciesMass(_biomassSp._eps[i]);
			}
			return m;
		}

		/**
		 * @return the volumetric fraction of eps in particle
		 */
		public float getEpsVolumetricFraction() {
			if (_biomassSp._eps == null)
				return 0;
			// sum volume occupied by all EPS components
			float epsVolume = 0;
			for (int i = 0; i < _biomassSp._eps.length; i++) {
				epsVolume += getSpeciesMass(_biomassSp._eps[i])
						/ _biomassSp._eps[i].getDensity();
			}
			// return the total volume occupied by eps by the total volume of
			// particle
			return epsVolume / _volume;
		}

		/**
		 * @return the value of relative growth rate (in relation to maximum
		 *         possible in system)
		 */
		public float getRelativeGrowth() {
			return _growthToMaximum;
		}

		/**
		 * @return the fraction of mass that corresponds to inert material
		 */
		public float getFractionOfInert() {
			float m = 0;
			float inert = 0;
			for (int i = 0; i < _masses.length; i++) {
				m += _masses[i];
				if (_particulateSp[i] != _biomassSp._inertMass) {
					inert = _masses[i];
				}
			}
			return inert / m;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Object#toString()
		 */
		public String toString() {
			String str = new String();
			for (int i = 0; i < _masses.length; i++) {
				str += _particulateSp[i].getName() + " - " + _masses[i] + "\n";
			}
			return str;
		}
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
		float[] masses = new float[_newParticleFractionalVolume.length];
		// the composition of each new cell
		// get the maximum volume of a particle
		float maxRadius = ((BiomassParticleContainer) Model.model().biomassContainer)
				.getMaximumRadius();
		// compute the maximum volume:
		float maxVolume = BiomassParticle.computeVolume(maxRadius);
		// compute the volume for this particle
		float volume = (0.8f + Model.model().getRandom() * 0.1f) * maxVolume;
		// set the mass valuesvalues in array
		for (int i = 0; i < masses.length; i++)
			masses[i] = volume * _newParticleFractionalVolume[i]
					* _particulates[i].getDensity();
		// call the instance specific inner class Composition
		// to allow more flexibility in extending BiomassSpecies
		return createComposition(this, masses);
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
		return new Composition(this, masses);
	}

	/**
	 * Return the array as a collection
	 * 
	 * @return collection of fixed substances
	 */
	public Collection<ParticulateSpecies> getFixedSpeciesAsCollection() {
		return Arrays.asList(_particulates);
	}

	/**
	 * @param s
	 *            The _activeMass to set.
	 * @throws ModelRuntimeException
	 *             if species is not contained in biomass species
	 */
	public void setActiveMass(ParticulateSpecies s)
			throws ModelRuntimeException {
		for (int i = 0; i < _particulates.length; i++) {
			if (s == _particulates[i]) {
				_activeMass = s;
				return;
			}
		}
		throw new ModelRuntimeException("Attempt to set non-existing species "
				+ s.getName() + " as active mass in " + getName());
	}

	/**
	 * Add a single species as EPS mass
	 * 
	 * @param mass
	 *            The _epsMass to set.
	 * @throws ModelException
	 *             if species is not contained in biomass species
	 */
	public void setEpsMass(ParticulateSpecies s) throws ModelException {
		ParticulateSpecies[] auxiliaryArray = { s };
		setEpsMass(auxiliaryArray);
	}

	/**
	 * Set an array of particulate species as the species composing the EPS
	 * capsule. Some rules: <bR>
	 * + All species in the s array must be already contained in the array of
	 * fixed species <br>
	 * + The s array must not contain repeated species. <br>
	 * + None of the species added must be already the active or inert mass
	 * species
	 * 
	 * @param s
	 *            the array of species to compose the EPS
	 * @throws ModelException
	 *             if any of the rules is broken
	 */
	public void setEpsMass(ParticulateSpecies[] s) throws ModelException {
		// cheack that no repetition exists
		try {
			assertThatNoRepetitionExists(s);
		} catch (ModelRuntimeException e) {
			throw new ModelException(e.toString());
		}
		// check that all species exist in _particulates array
		boolean[] checked = new boolean[s.length];
		for (int i = 0; i < s.length; i++) {
			// check is species is not active or inert already
			if ((s[i] == _activeMass) | (s[i] == _inertMass))
				throw new ModelException("Illegal attempt to" + " add species "
						+ s[i].getName() + " as part of eps mass in "
						+ getName()
						+ "; this species is already active or inert mass");
			// check that all species exist in _particulates array
			for (int j = 0; j < _particulates.length; j++) {
				if (s[i] == _particulates[j]) {
					checked[i] = true;
				}
			}
		}
		// check if all species passed the test where it is varified that they
		// exist on _particulates
		for (int i = 0; i < checked.length; i++) {
			if (!checked[i])
				throw new ModelException("Attempt to"
						+ " set non-existing species " + s[i].getName()
						+ " as part of eps mass in " + getName());
		}
		// if all checks where passed, the _eps array may be set
		_eps = s;
	}

	/**
	 * Make sure that no repetition of species exists in an array of particulate
	 * species
	 * 
	 * @param s
	 *            array to check
	 */
	private void assertThatNoRepetitionExists(ParticulateSpecies[] s) {
		// check if no species repetition exist
		for (int i = 0; i < s.length; i++) {
			for (int j = 0; j < s.length; j++) {
				if ((i != j) & (s[i] == s[j]))
					throw new ModelRuntimeException(
							"Attempt to add array of species"
									+ " with repeating species "
									+ s[i].getName() + " in " + getName());
			}
		}
	}

	/**
	 * @param mass
	 *            The _inertMass to set.
	 */
	public void setInertMass(ParticulateSpecies s) {
		for (int i = 0; i < _particulates.length; i++) {
			if (s == _particulates[i]) {
				_inertMass = s;
				return;
			}
		}
		throw new ModelRuntimeException("Attempt to set non-existing species "
				+ s.getName() + " as inert mass in " + getName());
	}

	/**
	 * @return the number of particles created in the present iteration
	 */
	public int getNumberOfParticlesCreatedInPresentIteration() {
		return _createdParticles;
	}

	/**
	 * Resets the counter for number of particles created in this iteration
	 */
	public void resetCreatedParticleCounter() {
		_createdParticles = 0;
	}

	/**
	 * Increments the counter for number of particles created in this iteration
	 */
	public void incrementCreatedParticleCounter() {
		++_createdParticles;
	}

	/**
	 * Individuals of this species will get their color from the relative growth
	 * at which they grow
	 */
	public void getColorFromGrowth() {
		_getColorFromGrowth = true;
	}

	/**
	 * true if biomass species has EPS species defined (i.e. _eps != null)
	 * 
	 * @return true if biomass species has eps species defined
	 */
	public boolean hasEpsCapsule() {
		return (_eps != null);
	}

	/**
	 * Check if a particulate species is contained in an array of particulate
	 * species
	 * 
	 * @param s
	 * @param array
	 * @return true is species s is contained in array, false otherwise
	 */
	private static boolean containedInArray(ParticulateSpecies s,
			ParticulateSpecies[] array) {
		if (array != null)
			for (int i = 0; i < array.length; i++) {
				if (s == array[i])
					return true;
			}
		return false;
	}

	/**
	 * Initialize the array of reactions. Get the reactions from the particulate
	 * species that compose this biomass species and store them in the array
	 * _reactions.
	 */
	private void initializeReactionsArray() {
		ArrayList allReactions = new ArrayList();
		// iterate through all the particulates and get the reactions
		// in which they are involved
		for (int i = 0; i < _particulates.length; i++) {
			List r;
			try {
				r = _particulates[i].getReactionsInvlovedAsArrayList();
			} catch (NullPointerException e) {
				// TODO Auto-generated catch block
				throw new ModelRuntimeException(
						"Check if you associated processes with species "
								+ _particulates[i].getName());
			}
			for (Iterator iter = r.iterator(); iter.hasNext();) {
				Reaction element = (Reaction) iter.next();
				// ensure that reaction was not added already
				if (!allReactions.contains(element))
					allReactions.add(element);
			}
		}
		// initialize the _reactions array and fill it
		_reactions = new Reaction[allReactions.size()];
		for (int i = 0; i < _reactions.length; i++) {
			_reactions[i] = (Reaction) (allReactions.get(i));
		}
	}

	/**
	 * Get the shoving hierarchy. Species with higher hierarchy push ones with
	 * lower, but not vice-versa
	 * 
	 * @return the shoving hierarchy
	 */
	public int getShovingHierarchy() {
		return _shovingHierarchy;
	}

	/**
	 * Get the shoving hierarchy. Species with higher hierarchy push ones with
	 * lower, but not vice-versa
	 */
	public void setShovingHierarchy(int shovingHierarchy) {
		_shovingHierarchy = shovingHierarchy;
	}

	@Override
	public String toString() {
		return  _name;
	}
	
	
}