package nl.tudelft.bt.model.work.epsproducers.sinkinganalogy;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

/**
 * Implements the mutator particle
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Oct 18, 2006
 */
public class MutatorBiomassParticle extends BiomassParticle {
	private MutatorBiomassSpecies _species;

	public static float MUTATIONPROB = 0.1f;

	private static boolean isWildType;

	public MutatorBiomassParticle(MutatorBiomassSpecies s) {
		super(s);
		_species = s;
		// start simulation as wild-type
		isWildType = true;
	}

	@Override
	public BiomassParticle divide() {
		MutatorBiomassParticle daughter = (MutatorBiomassParticle) super
				.divide();
		// code for mutation
		// this ensures that mutation only occurs once per simulation
		if (Model.model().getRandom() < MUTATIONPROB) {
			mutate();
		}
		return daughter;
	}

	@Override
	protected Object clone() throws CloneNotSupportedException {
		return super.clone();
	}

	/**
	 * execute the mutation if randonmly chosen for this
	 */
	private void mutate() {
		isWildType = !isWildType;
		if (isWildType) {
			// pass all mass from the wild-type trait to the mutant trait
			_composition.setMass(_species._mutantTrait, _composition
					.getSpeciesMass(_species._wildTypeTrait));
			// and set the WT to 0
			_composition.setMass(_species._wildTypeTrait, 0);
		} else {
			// pass all mass from the mutant traite trait to wild-type
			_composition.setMass(_species._wildTypeTrait, _composition
					.getSpeciesMass(_species._mutantTrait));
			// and set the WT to 0
			_composition.setMass(_species._mutantTrait, 0);
		}
	}
}
