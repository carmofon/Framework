package nl.tudelft.bt.model.work.epsproducers.mutation;

import java.awt.Color;
import java.util.Collection;
import java.util.Iterator;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

/**
 * Implements the mutator particle
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Oct 18, 2006
 */
public class MutatorBiomassParticle extends BiomassParticle {
	private MutatorBiomassSpeciesForImpact _species;

	private static boolean _mutated = false;

	public MutatorBiomassParticle(MutatorBiomassSpeciesForImpact s) {
		super(s);
		_species = s;
	}

	@Override
	public BiomassParticle divide() {
		MutatorBiomassParticle daughter = (MutatorBiomassParticle) super
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
					// ///// change all the neighbors
					// iterate thgouh collection of particles
					Collection particles = _particleContainer
							.getBiomassAsBiomassParticleCollection();
					for (Iterator i = particles.iterator(); i.hasNext();) {
						BiomassParticle p = (BiomassParticle) i.next();
						// check if particle belongs to this species
						if (p != this) {
							if (p.isOfSpecies(_species)) {
								MutatorBiomassParticle p2 = ((MutatorBiomassParticle) p);
								if (this.distanceTo(p2) <= _species._neighborhoodRadius) {
									p2.changeToNeighbor();
								}
							}
						}
					}
					daughter.changeToNeighbor();
				}
		return daughter;
	}

	@Override
	protected Object clone() throws CloneNotSupportedException {
		// TODO Auto-generated method stub
		return super.clone();
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
	}

	/**
	 * execute the change to neighbor trait
	 */
	private void changeToNeighbor() {
		// pass all mass from the wild-type trait to the mutant trait
		_composition.setMass(_species._neighborsTrait, _composition
				.getSpeciesMass(_species._wildTypeTrait));
		// and set the WT to 0
		_composition.setMass(_species._wildTypeTrait, 0);
	}
}
