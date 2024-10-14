const express = require("express");
const mongoose = require("mongoose");
const rasInhibitor = require("./models/rasInhibitor");
const cors = require("cors");
const bodyParser = require("body-parser");

const main = async () => {
  await mongoose.connect("mongodb://localhost:27017/grad_school");
  console.log("connected to mongo");
};

const wrapper = () => {
  const app = express();
  app.use(cors());
  app.use(bodyParser.json());
  app.get("/", (req, res) => {
    console.log(req.body);
    console.log(req);
    console.log(req.headers);
    res.sendStatus(200);
  });
  app.get("/ras-inhibitors", async (req, res) => {
    res.json(await rasInhibitor.find({}));
  });
  app.post("/remove-pdb", async (req, res) => {
    try {
      const { drugId, pdbId } = req.body;
      console.log(drugId, pdbId);
      if (!drugId || !pdbId) {
        console.log("drugId or pdbId not provided");
        return res.status(400).send("drugId or pdbId not provided");
      }
      const result = await rasInhibitor.findById(drugId);
      result.pdbHits = result.pdbHits.filter((hit) => hit != pdbId);
      await rasInhibitor.findByIdAndUpdate(drugId, { pdbHits: result.pdbHits });
      res.sendStatus(200);
    } catch (err) {
      console.log(err);
      res.status(400).send(err);
    }
  });
  app.post("/add-note", async (req, res) => {
    try {
      const { drugId, newNote } = req.body;
      console.log(drugId, newNote);
      if (!drugId || !newNote) {
        console.log("drugId or newNote were empty");
        return res.status(400).send("drugId or newNote were empty");
      }
      const result = await rasInhibitor.findById(drugId);
      result.notes.push(newNote);
      await rasInhibitor.findByIdAndUpdate(drugId, { notes: result.notes });
      res.sendStatus(200);
    } catch (err) {
      res.status(400).send(err);
    }
  });
  app.post("/delete-note", async (req, res) => {
    try {
      const { drugId, oldNote } = req.body;
      console.log(drugId, oldNote);
      if (!drugId || !oldNote) {
        console.log("drugId or oldNote not provided");
        return res.status(400).send("drugId or oldNote not provided");
      }
      const result = await rasInhibitor.findById(drugId);
      result.notes = result.notes.filter((note) => note != oldNote);
      await rasInhibitor.findByIdAndUpdate(drugId, { notes: result.notes });
      res.sendStatus(200);
    } catch (err) {
      console.log(err);
      res.status(400).send(err);
    }
  });
  app.listen("8080", () => {
    console.log("listening on 8080");
  });
};

main().then(wrapper);
