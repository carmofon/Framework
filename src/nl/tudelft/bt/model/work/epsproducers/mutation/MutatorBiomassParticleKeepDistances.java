package nl.tudelft.bt.model.work.epsproducers.mutation;

import java.awt.Color;
import java.util.Collection;
import java.util.Iterator;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.util.ColorMaps;
import nl.tudelft.bt.model.util.ExtraMath;

/**
 * Implements the mutator particle. A mutation will occur at a given time. The
 * change in trait will be defined by the biomass species. The mutante will
 * change its trait and will tag the remaining particles with the distance at
 * which they are to the mutant.
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Oct 18, 2006
 */
public class MutatorBiomassParticleKeepDistances extends BiomassParticle {
	private BiomassSpeciesForKeepingDistances _species;

	// tag that allows to quickly determine if mutation has already occurred
	private static boolean _mutated = false;

	//
	private boolean _isMutant = false;

	private float _distanceToMutant;

	private static float _maximumDistanceToMutant;

	public MutatorBiomassParticleKeepDistances(
			BiomassSpeciesForKeepingDistances s) {
		super(s);
		_species = s;
	}

	@Override
	public BiomassParticle divide() {
		MutatorBiomassParticleKeepDistances daughter = (MutatorBiomassParticleKeepDistances) super
				.divide();
		// code for mutation
		// this ensures that mutation only occurs once per simulation
		if (!isEpsOnly())
			if (!_mutated)
				// determine if particle must mutate
				if (Model.model().getTime() > _species._mutationTime) {
					mutate();
					//
					_mutated = true;
					_maximumDistanceToMutant = 0;
					// ///// change all the neighbors
					// iterate thgouh collection of particles and tag them with
					// distance to this particle
					Collection particles = _particleContainer
							.getBiomassAsBiomassParticleCollection();
					for (Iterator i = particles.iterator(); i.hasNext();) {
						BiomassParticle p = (BiomassParticle) i.next();

						if (p != this) {
							// check if particle belongs to this species
							if (p.isOfSpecies(_species)) {
								MutatorBiomassParticleKeepDistances p2 = ((MutatorBiomassParticleKeepDistances) p);
								if (!p2.isEpsOnly()) {
									p2._distanceToMutant = distanceTo(p2);
									// update the maximum distance to mutant
									_maximumDistanceToMutant = ExtraMath.max(
											_maximumDistanceToMutant,
											p2._distanceToMutant);
								}
							}
						}
					}
					daughter._distanceToMutant = distanceTo(daughter);
				}
		return daughter;
	}

	public boolean isEpsOnly() {
		return _composition.getTotalMass() == _composition.getEpsMass();
	}

	/**
	 * execute the mutation if randonmly chosen for this
	 */
	private void mutate() {
		// pass all mass from the wild-type trait to the mutant trait
		_composition.setMass(_species._mutantTrait, _composition
				.getSpeciesMass(_species._wildTypeTrait));
		// and set the WT to 0
		_composition.setMass(_species._wildTypeTrait, 0);
		// and the flag to mutant
		_isMutant = true;
		_distanceToMutant = 0;
	}

	/**
	 * @return true if this is the mutant
	 */
	public boolean isMutant() {
		return _isMutant;
	}

	/**
	 * @return the distance to mutant of ancestor cell at the time of mutation
	 */
	public float getDistanceToMutant() {
		return _distanceToMutant;
	}

	@Override
	public Color getColorCore() {
		if ((_mutated) & (!_isMutant))
			return ColorMaps.getJetColor(_distanceToMutant
					/ _maximumDistanceToMutant);
		return super.getColorCore();
	}

	@Override
	protected Object clone() throws CloneNotSupportedException {
		// TODO Auto-generated method stub
		return super.clone();
	}

	/**
	 * @return true is mutation already occurred
	 */
	public static boolean mutationOccurred() {
		return _mutated;
	}
	
	/**
	 * @return the total mass minus the EPS mass
	 */
	public float getActiveMass() {
		return getTotalMass() - getTotalEps();
	}
}
