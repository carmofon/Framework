package nl.tudelft.bt.model.work.granule;

import nl.tudelft.bt.model.bulkconcentrations.BulkConcentration;
import nl.tudelft.bt.model.exceptions.ModelException;

public class AnaerobicFeedControlledDO40Precision0_01MaxF0_95_eta05 extends UpdatedAquasimModel {

	protected void defineSpeciesAndReactions() throws ModelException {
		precision = 0.01f;
		maxFractionToDecrease = 0.95f;
		super.defineSpeciesAndReactions();
	}

	public BulkConcentration createOxygen() {
		oxygenBulkConcentration = 8e-3f * 0.4f;
		float aerationStart = 1.0f;
		float tresholdValue = ammoniumFeedConcentration * 0.01f;
		return new ItermittentBulkConcentrationControlled(oxygenBulkConcentration,
				aerationStart, _cycle,ammonium, nitrite,
				tresholdValue);
	}

	public static void main(String[] args) {
		AnaerobicFeedControlledDO40Precision0_01MaxF0_95_eta05 app = new AnaerobicFeedControlledDO40Precision0_01MaxF0_95_eta05();
		app.run(args);
	}
}
