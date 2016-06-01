package nl.tudelft.bt.model.work.will;

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
	private MutatorBiomassSpecies _species;

	private boolean _mutated = false;

	public MutatorBiomassParticle(MutatorBiomassSpecies s) {
		super(s);
		_species = s;
	}

	@Override
	public BiomassParticle divide() {
		MutatorBiomassParticle daughter = (MutatorBiomassParticle) super
				.divide();
		// code for mutation
		// this ensures that mutation only occurs once per simulation
		if (!_mutated)
			// determine if particle must mutate
			if (_species._mutationRate > Model.model().getRandom()) {
				mutate();
				//
				_mutated = true;
			
			}
		return daughter;
	}

	@Override
	protected Object clone() throws CloneNotSupportedException {
		return super.clone();
	}

	/**
	 * execute the mutation if randomly chosen for this
	 */
	private void mutate() {
		// pass all mass from the wild-type trait to the mutant trait
		_composition.setMass(_species._mutantTrait, _composition
				.getSpeciesMass(_species._wildTypeTrait));
		// and set the WT to 0
		_composition.setMass(_species._wildTypeTrait, 0);
	}
}
