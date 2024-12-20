// const mysql = require("mysql2");
// const fs = require("fs");
// const path = require("path");
const initRDKitModule = require("@rdkit/rdkit");

// const delay = (delayInms) => {
//   return new Promise((resolve) => setTimeout(resolve, delayInms));
// };

// const createInchiKeys = async () => {
//   const allMonomers = [];
//   const monomersToUpdate = [];
//   let shouldContinue = true;
//   const monomerQuery = (offset, limit) =>
//     new Promise(async (resolve, reject) => {
//       const rdkit = await initRDKitModule();
//       const connection = mysql.createConnection({
//         uri: "jdbc:mysql://localhost:3306",
//         host: "localhost",
//         user: "root",
//         database: "binding_data",
//       });
//       connection.query(
//         `SELECT x.* FROM binding_data.monomer x LIMIT ${offset},${limit}`,
//         (err, res, fields) => {
//           if (err) {
//             shouldContinue = false;
//             throw new Error(err);
//             return;
//           }
//           if (!res?.length) {
//             shouldContinue = false;
//             return;
//           }
//           for (result of res) {
//             if (result.inchi_key && result.inchi) continue;
//             let mol;
//             let smiles;
//             try {
//               smiles = result.smiles_string;
//               mol = rdkit.get_mol(smiles);
//             } catch {
//               console.log("could not construct molecule");
//               continue;
//             }
//             let inchi = "";
//             let inchi_key = "";
//             try {
//               inchi = mol.get_inchi();
//               inchi_key = rdkit.get_inchikey_for_inchi(inchi);
//             } catch {
//               inchi = "INCHI NOT FOUND";
//               inchi_key = "INCHI KEY NOT FOUND";
//             }
//             monomersToUpdate.push({ ...result, inchi, inchi_key });
//           }
//           const prevData = JSON.parse(
//             fs.readFileSync(path.join(__dirname, "tempData.json"))
//           );
//           prevData.push(...monomersToUpdate);
//           fs.writeFileSync(
//             path.join(__dirname, "tempData.json"),
//             JSON.stringify(monomersToUpdate)
//           );
//           resolve();
//           connection.end();
//           delete rdkit;
//           return;
//         }
//       );
//     });
//   let i = 0;
//   while (shouldContinue) {
//     console.log("round #" + i);
//     query = monomerQuery(i * 1000, (i + 1) * 1000);
//     await query;

//     i += 1;
//   }
// };

// createInchiKeys();

const mysql = require("mysql2/promise");

const processSmilesString = async (smiles, rdkit) => {
  let mol;
  let inchi;
  let inchiKey;
  let shouldRestartRdkit = false;
  try {
    mol = rdkit.get_mol(smiles);
    inchi = mol.get_inchi(); // Example output
    inchiKey = rdkit.get_inchikey_for_inchi(inchi); // Example output
  } catch {
    try {
      console.log("restarting rdkit");
      shouldRestartRdkit = true;
      rdkit = await initRDKitModule();
      mol = rdkit.get_mol(smiles);
      inchi = mol.get_inchi(); // Example output
      inchiKey = rdkit.get_inchikey_for_inchi(inchi); // Example output
    } catch {
      inchi = "COULD NOT CALCULATE INCHI";
      inchiKey = "COULD NOT FIND INCHI_KEY";
    }
    inchi = "COULD NOT CALCULATE INCHI";
    inchiKey = "COULD NOT FIND INCHI_KEY";
  }
  // Dummy function to process smiles_string.
  // Replace this with your actual processing logic.

  return { inchi, inchiKey, shouldRestartRdkit };
};

const updateEntries = async (connection, entries) => {
  let rdkit = await initRDKitModule();
  for (const entry of entries) {
    const { monomerid, smiles_string } = entry;
    if (!entry || !smiles_string) continue;
    console.log(monomerid, smiles_string);
    const { inchi, inchiKey, shouldRestartRdkit } = await processSmilesString(
      smiles_string,
      rdkit
    );
    if (shouldRestartRdkit) rdkit = await initRDKitModule();

    // Update the database
    await connection.execute(
      "UPDATE monomer SET inchi = ?, inchi_key = ? WHERE monomerid = ?",
      [inchi, inchiKey, monomerid]
    );
  }
};

const main = async () => {
  const connection = await mysql.createConnection({
    host: "localhost",
    user: "root", // replace with your MySQL username
    database: "binding_data",
  });

  const batchSize = 1000; // Number of records to process in each batch
  let offset = 645336;

  try {
    while (true) {
      const [rows] = await connection.execute(
        `SELECT * FROM binding_data.monomer LIMIT ${offset},${batchSize};`
      );

      if (rows.length === 0) {
        break; // Exit loop if no more records
      }
      await updateEntries(connection, rows);
      offset += batchSize;

      console.log(`Processed ${offset} records...`);
    }
  } catch (error) {
    console.error("Error processing records:", error);
  } finally {
    await connection.end();
  }
};

main().catch(console.error);
