const mongoose = require("mongoose");

module.exports = mongoose.model(
  "ras_inhibitors",
  new mongoose.Schema({
    name: String,
    url: String,
    prices: [{ amount: String, price: Number }],
    image_url: String,
    ccdc_url: String,
    pdb_url: String,
    notes: [String],
    smiles: String,
    description: String,
    synonyms: [String],
    hits: [{ name: String, hits: [String] }],
    pdbHits: [String],
  })
);
