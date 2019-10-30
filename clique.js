/*
 *  Clique Generation with Node.Js
 */
const fs = require('fs');
const readline = require('readline');

const inputReader = readline.createInterface({
  input: process.stdin,
  output: process.stdout
});

console.log('Enter the order of the clique,\
 or enter \'q\' to quit the program')

inputReader.on('line', (line) => {
  let n = Number(line);
  if (!Number.isInteger(n))
    process.exit();

  let lines = '';
  lines += `${n}\n`;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      lines += j === i ? '' : `${j} `;
    }
    lines += '\n';
  }

  fs.writeFile("./alists/clique.alist", lines, function (err) {
    if (err) {
      return console.log(err);
    }
    console.log("The file was saved!");
  });
});