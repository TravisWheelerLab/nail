import * as soda from "https://esm.run/@sodaviz/soda@0.13.1";

function run() {
  let chart = new soda.Chart({
    selector: "div.container",
    zoomable: true,
    resizable: true,
    rowHeight: 30,
    draw(params) {
      this.addAxis();

      soda.chevronRectangle({
        chart: this,
        annotations: params.forward.leftBounds,
        row: (d) => d.a.row,
        orientation: soda.Orientation.Forward,
        strokeColor: "blue",
        strokeOpacity: 0.5,
      })

      soda.chevronRectangle({
        chart: this,
        annotations: params.forward.rightBounds,
        row: (d) => d.a.row,
        orientation: soda.Orientation.Forward,
        strokeColor: "red",
        strokeOpacity: 0.5,
      })

      soda.chevronRectangle({
        chart: this,
        annotations: params.backward.leftBounds,
        row: (d) => d.a.row,
        orientation: soda.Orientation.Reverse,
        strokeColor: "blue",
        strokeOpacity: 0.5,
      })

      soda.chevronRectangle({
        chart: this,
        annotations: params.backward.rightBounds,
        row: (d) => d.a.row,
        orientation: soda.Orientation.Reverse,
        strokeColor: "red",
        strokeOpacity: 0.5,
      })

      soda.line({
        chart: this,
        annotations: params.forward.antiDiagonals,
        x1: (d) => d.c.xScale(d.a.x1),
        x2: (d) => d.c.xScale(d.a.x2),
        y1: (d) => d.c.rowHeight * d.a.y1,
        y2: (d) => d.c.rowHeight * d.a.y2,
        strokeWidth: 2,
        strokeColor: "green",
        strokeOpacity: 0.5,
      });

      soda.line({
        chart: this,
        annotations: params.backward.antiDiagonals,
        x1: (d) => d.c.xScale(d.a.x1),
        x2: (d) => d.c.xScale(d.a.x2),
        y1: (d) => d.c.rowHeight * d.a.y1,
        y2: (d) => d.c.rowHeight * d.a.y2,
        strokeWidth: 2,
        strokeColor: "red",
        strokeOpacity: 0.5,
      });

    }
  });

  function prepareBounds(boundStrings, prefix) {
    let left = boundStrings[0].split("|");
    let right = boundStrings[1].split("|");

    let rowColFn = (s) => {
      let tokens = s.split(",");
      let row = parseInt(tokens[0]);
      let col = parseInt(tokens[1]);
      return [row, col];
    };

    let leftBounds = [];
    let rightBounds = [];
    let antiDiagonals = [];
    for (let i = 0; i < left.length; i++) {
      let [leftRow, leftCol] = rowColFn(left[i]);
      let [rightRow, rightCol] = rowColFn(right[i]);
      let antiDiagonal = leftRow + leftCol;

      leftBounds.push({
        id: `${prefix}l-${antiDiagonal}`,
        antiDiagonal,
        row: leftRow,
        col: leftCol,
        start: leftCol,
        end: leftCol + 1,
      })

      rightBounds.push({
        id: `${prefix}r-${antiDiagonal}`,
        antiDiagonal,
        row: rightRow,
        col: rightCol,
        start: rightCol,
        end: rightCol + 1,
      })

      antiDiagonals.push({
        id: `${prefix}a-${antiDiagonal}`,
        antiDiagonal,
        x1: leftCol + 0.5,
        x2: rightCol + 0.5,
        y1: leftRow + 0.5,
        y2: rightRow + 0.5,
      })

    }

    return { leftBounds, rightBounds, antiDiagonals };
  }

  let forward =  prepareBounds(data.forwardBounds, "f");
  let backward = prepareBounds(data.backwardBounds, "b");
  let end = 0;
  let rowCount = 0;
  
  for (const f of forward.rightBounds.concat(forward.leftBounds)) {
    end = Math.max(end, f.end + 1);
    rowCount = Math.max(rowCount, f.row + 1)
  }

  for (const b of backward.leftBounds.concat(backward.rightBounds)) {
    end = Math.max(end, b.end + 1);
    rowCount = Math.max(rowCount, b.row + 1)
  }
  
  let params = {
    start: 0,
    end,
    rowCount,
    forward,
    backward,
  };
  console.log(params); 

  chart.render(params);
}
