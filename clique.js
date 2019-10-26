let n = 25;
console.log(n);
for (let i = 0; i < n; i++) {
  let line = '';
  for (let j = 0; j < n; j++) {
    line += j === i ? '' : `${j} `;
  }
  console.log(line);
}