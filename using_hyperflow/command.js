var spawn = require('child_process').spawn;

function command(ins, outs, executor, config, cb) {
    var exec = config.executor.executable,
        args = config.executor.args;

    console.log("Executing:", exec, args);

    var proc = spawn(exec, [ args ]);

    proc.stdout.on('data', function(data) {
        console.log(exec, 'stdout:' + data);
    });

    proc.stderr.on('data', function(data) {
        console.log(exec, 'stderr:' + data);
    });

    proc.on('exit', function(code) {
        console.log(exec, 'exiting with code:' + code);
        cb(null, outs);
    });

    proc.on('close', function (code, signal) {
        console.log(exec, 'terminated due to receipt of signal '+signal);
    });
}

exports.command = command;

// copied from http://stackoverflow.com/questions/105034/how-to-create-a-guid-uuid-in-javascript
function createUUID() {
    // http://www.ietf.org/rfc/rfc4122.txt
    var s = new Array(36);
    var hexDigits = "0123456789abcdef";
    for (var i = 0; i < 36; i++) {
        s[i] = hexDigits.substr(Math.floor(Math.random() * 0x10), 1);
    }
    s[14] = "4";  // bits 12-15 of the time_hi_and_version field to 0010
    s[19] = hexDigits.substr((s[19] & 0x3) | 0x8, 1);  // bits 6-7 of the clock_seq_hi_and_reserved to 01
    s[8] = s[13] = s[18] = s[23] = "-";

    var uuid = s.join("");
    return uuid;
}

