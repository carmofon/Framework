package nl.tudelft.bt.model.work.granule;

import nl.tudelft.bt.model.bulkconcentrations.BulkConcentration;
import nl.tudelft.bt.model.bulkconcentrations.ConstantBulkConcentration;
import nl.tudelft.bt.model.exceptions.ModelException;

public class FullAerobicDO40Precision0_01MaxF0_95_eta05 extends UpdatedAquasimModel {

	protected void defineSpeciesAndReactions() throws ModelException {
		precision = 0.01f;
		maxFractionToDecrease = 0.95f;
		super.defineSpeciesAndReactions();
	}

	public BulkConcentration createOxygen() {
		oxygenBulkConcentration = 8e-3f * 0.4f;
		return new ConstantBulkConcentration(oxygenBulkConcentration);
	}

	public static void main(String[] args) {
		FullAerobicDO40Precision0_01MaxF0_95_eta05 app = new FullAerobicDO40Precision0_01MaxF0_95_eta05();
		app.run(args);
	}
}
