package nl.tudelft.bt.model.work.vanni;

import java.awt.Color;
import java.util.Collection;
import java.util.Iterator;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

/**
 * Implements the mutator particle
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Oct 18, 2006
 */
public class MutatorBiomassParticle extends BiomassParticle {
	private MutatorBiomassSpecies _species;

	private int mutationLevel= 0;
	
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
			// determine if particle must mutate
			if (_species._mutationRate > Model.model().getRandom()) {
				mutate();
			
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
		if (mutationLevel == _species._numberOfLevels-1) return;
		//Do mutation
		ParticulateSpecies presentLevel = _species.arraySpecies[mutationLevel];
		ParticulateSpecies nextLevel = _species.arraySpecies[mutationLevel+1];
		mutationLevel++;
		// pass all mass from the wild-type trait to the mutant trait
		_composition.setMass(nextLevel, _composition
				.getSpeciesMass(presentLevel));
		// and set the WT to 0
		_composition.setMass(presentLevel, 0);
	}
	
	
	public Color getColorCore() {
		return super.getColorCore();
		/*if (this.mutationLevel > 0) {
			int trash = 0;
			trash++;
		}
		Color c = _species.arraySpecies[this.mutationLevel].getColor();
		return c;*/
	}
}
