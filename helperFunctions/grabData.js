const { parse } = require("node-html-parser");
const mongoose = require("mongoose");
const rasInhibitor = require("../models/rasInhibitor");
const fs = require("fs");
const path = require("path");
const grabNamesMedChemExpress = async () => {
  await mongoose.connect("mongodb://localhost:27017/grad_school");

  mongoose.model("ras_inhibitors", rasInhibitor);
  const rasInhibitors = mongoose.model("ras_inhibitors");
  const drugNames = await fetch(
    "https://www.medchemexpress.com/Targets/Ras/Ras-comparison.html"
  );
  const webPage = parse(await drugNames.text());
  const tableData = webPage.querySelectorAll(".col-pro-name");
  for (el of tableData) {
    if (!el.querySelector("a")) continue;
    const { innerHTML, _attrs } = el.querySelector("a");
    await rasInhibitors.create({ name: innerHTML, url: _attrs.href });
  }
  mongoose.connection.close();
};

const grabDrugDataMedChemExpress = async () => {
  await mongoose.connect("mongodb://localhost:27017/grad_school");
  mongoose.model("ras_inhibitors", rasInhibitor);
  const rasInhibitors = mongoose.model("ras_inhibitors");

  const drugs = await rasInhibitors.find({});

  for (let i = 0; i < drugs.length; i++) {
    const drug = drugs[i].toObject();

    const data = await fetch("https://www.medchemexpress.com" + drug.url);

    const webPage = parse(await data.text());
    const masses = webPage.querySelectorAll(".pro_price_1");
    const prices = webPage.querySelectorAll(".pro_price_2");
    const image_url = webPage.querySelector(".data-img");
    const smiles = webPage.querySelectorAll("tr");
    const description = webPage.querySelector("#product_syn");

    const cleanedPrices = [];
    for (let i = 0; i < masses.length; i++) {
      if (masses[i]?.innerHTML.includes("Size")) continue;
      const mass = masses[i]?.innerHTML?.replace(/<[^<|>]*>/g, "").trim();
      let price = parseInt(
        prices[i]?.querySelector("span")?.innerHTML?.replace("USD", "").trim()
      );
      if (isNaN(price)) price = null;
      cleanedPrices.push({ amount: mass, price });
    }

    const cleanedSmiles = smiles
      .filter((el) => el.innerHTML.includes("SMILES"))[0]
      ?.querySelector("p")?.innerHTML;
    const cleanedImageUrl = image_url?._attrs.src.replace("//", "");
    const cleanedDescription = description?.innerHTML
      .replace(/<[^<|>]*>/g, "")
      .trim();
    console.log(
      cleanedPrices,
      cleanedImageUrl,
      cleanedSmiles,
      cleanedDescription
    );
    const result = await rasInhibitors.findOneAndUpdate(
      { name: drug.name },
      {
        $set: {
          image_url: cleanedImageUrl,
          prices: cleanedPrices,
          smiles: cleanedSmiles,
          description: cleanedDescription,
        },
      },
      { new: true }
    );
    await delay(1000);
  }
  mongoose.connection.close();
};

const grabSynonyms = async () => {
  await mongoose.connect("mongodb://localhost:27017/grad_school");
  const results = await rasInhibitor.find({});

  for (let i = 0; i < results.length; i++) {
    const result = results[i].toObject();
    console.log("getting synonyms for " + result?.name + "with SMILES");
    let synonyms;
    let synonymsData;
    try {
      synonyms = await fetch(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/" +
          encodeURI(result.smiles) +
          "/synonyms/JSON"
      );
      synonymsData = await synonyms.json();
    } catch {
      console.log("SMILES search failed, trying with compound name");
      try {
        await delay(1000);
        synonyms = await fetch(
          "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" +
            encodeURI(result.name) +
            "/synonyms/JSON"
        );
        synonymsData = await synonyms.json();
      } catch {
        console.log("name search failed, defaulting to empty array");
        await rasInhibitor.findOneAndUpdate(
          { name: result.name },
          { synonyms: [] }
        );
        await delay(1000);
        continue;
      }
    }

    await rasInhibitor.findOneAndUpdate(
      { name: result.name },
      { synonyms: synonymsData.InformationList?.Information[0]?.Synonym }
    );
    await delay(1000);
  }
  console.log("done");
};

const delay = (delayInms) => {
  return new Promise((resolve) => setTimeout(resolve, delayInms));
};

const collateSynonyms = async () => {
  await mongoose.connect("mongodb://localhost:27017/grad_school");
  const rasInhibs = await rasInhibitor.find({});
  const collatedSynonyms = [];
  for (let i = 0; i < rasInhibs.length; i++) {
    rasInhib = rasInhibs[i].toObject();
    console.log(rasInhib.synonyms);
    collatedSynonyms.push({
      drugName: rasInhib.name,
      synonyms: rasInhib.synonyms,
      smiles: rasInhib.smiles,
    });
  }
  fs.writeFileSync(
    path.join(__dirname, "../data/collatedSynonyms.json"),
    JSON.stringify(collatedSynonyms)
  );
};

const addCSDHitsToDB = async () => {
  await mongoose.connect("mongodb://localhost:27017/grad_school");
  let data = fs.readFileSync(
    path.join(__dirname, "../data/collatedSynonyms_withHits.json")
  );
  data = JSON.parse(data);
  for (let i = 0; i < data.length; i++) {
    const datum = data[i];
    const hits = [];
    for (let j = 0; j < datum?.hits?.length; j++) {
      const [name, hit] = Object.entries(datum.hits[j])[0];
      hits.push({ name, hits: hit });
    }
    await rasInhibitor.findOneAndUpdate({ name: datum.drugName }, { hits });
  }
  await mongoose.connection.close();
};

const addPDBHitsToDB = async () => {
  await mongoose.connect("mongodb://localhost:27017/grad_school");
  let data = fs.readFileSync(
    path.join(__dirname, "../data/collatedSynonyms_withPDBHits.json")
  );
  data = JSON.parse(data);
  for (let i = 0; i < data.length; i++) {
    const { pdbHits, drugName } = data[i];
    await rasInhibitor.findOneAndUpdate({ name: drugName }, { pdbHits });
  }
  await mongoose.connection.close();
};

const dumpData = async () => {
  await mongoose.connect("mongodb://localhost:27017/grad_school");
  const data = await rasInhibitor.find({});
  fs.writeFileSync(
    path.join(__dirname, "../data/db_dumpjson"),
    JSON.stringify(data)
  );
  return;
};

dumpData();
