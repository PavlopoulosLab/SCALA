// This function creates an alert message for the user
// @param messsage: Error message to be printed
// @return true
function shinyAlert(message){
  alert(message);
  return true;
}

// This function sends a string to the browser console log
// @param messsage: Message to be printed
// @return true
function shinyLog(message){
  console.log(message);
  return true;
}

// This function disappears the requested button so the user does not
// click it multiple times before current function execution is finished
// @param btn: button element id
function disableButton(btn){
  var button = document.getElementById(btn);
  button.style.display = "none";
  return true;
}

// This function re-appears the requested button
// after the end of fucntion execution
// @param btn: button element id
function enableButton(btn){
  var button = document.getElementById(btn);
  button.style.display = "inline-block";
  return true;
}

// This function starts a requested loader
// @param m[0]: loader id
// @param m[1]: percentage to fill loader
function startLoader(m){
  var loader = document.getElementById(m[0]);
  loader.nextElementSibling.style.opacity = 0.5;
  loader.style.opacity = 1;
  loader.style.display = "inline-block";
  var bar = loader.ldBar;
  bar.set(parseInt(m[1]));
  return true;
}

// This function ends a requested loader
// @param m: loader id
function finishLoader(m){
  var loader = document.getElementById(m);
  var bar = loader.ldBar;
  bar.set(0);
  loader.nextElementSibling.style.opacity = 1;
  loader.style.display = "none";
  return true;
}

// This function calculates the height of a given div
// @param m[0]: inner div id
// @param m[1]: #plots for Manhatan
function fixHeight(m){
  var height_unit = 400,
      div = document.getElementById(m[0]).parentElement.parentElement.parentElement.parentElement.parentElement;
  div.style.height = String(height_unit * m[1]).concat("px");
  return true;
}

Shiny.addCustomMessageHandler("handler_alert", shinyAlert);
Shiny.addCustomMessageHandler("handler_log", shinyLog);
Shiny.addCustomMessageHandler("handler_disableButton", disableButton);
Shiny.addCustomMessageHandler("handler_enableButton", enableButton);
Shiny.addCustomMessageHandler("handler_startLoader", startLoader);
Shiny.addCustomMessageHandler("handler_finishLoader", finishLoader);
Shiny.addCustomMessageHandler("handler_fixHeight", fixHeight);
