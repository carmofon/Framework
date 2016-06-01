package nl.tudelft.bt.model.work.will;

import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

public class AdhesiveBiomassSpecies extends ChemotacticBiomassSpecies {

	private float _ws;
	
	private float _wg;
	
	private float _ra;
	
	/**
	 * @param name
	 * @param species
	 * @param fractionalCompositionInVolume
	 * @param attractant
	 * 					the chemoattractant solute species
	 * @param motilityDiffusivity
	 * @param Kspeed
	 * 					migration speed half-max constant [g/L]
	 * @param Kspread
	 * 					migration directedness half-max constant [um^-1]
	 * @param ws
	 * 					weight coefficient of social factors on migration [dimless]
	 * @param wg
	 * 					weight coefficient of individual factors on migration [dimless]
	 * @param relativeAdhesiveRadius
	 * 					effective adhesion range of cell relative to cell radius
	 * @throws NonMatchingNumberException
	 */
	
	public AdhesiveBiomassSpecies(String name, ParticulateSpecies[] species,
			float[] fractionalCompositionInVolume, SoluteSpecies attractant,
			float motilityDiffusivity, float Kspeed, float Kspread, float minMotility,
			float ws, float wg, float relativeAdhesiveRadius)
			throws NonMatchingNumberException {
		super(name, species, fractionalCompositionInVolume, attractant,
				motilityDiffusivity, Kspeed, Kspread, minMotility);
		// require that the social and individual weights be normalized
//		if (ws + wg != 1) {
//			throw new NonMatchingNumberException(
//					"Social and individual weights must sum to 1");
//		}
		_ws = ws;
		_wg = wg;
		_ra = relativeAdhesiveRadius;
	}
	
	public float getws() {
		return _ws;
	}
	
	public float getwg() {
		return _wg;
	}
	
	public float getra() {
		return _ra;
	}
	
	public BiomassParticle createBiomassParticle() {
		return new AdhesiveBiomassParticle(this);
	}

}
