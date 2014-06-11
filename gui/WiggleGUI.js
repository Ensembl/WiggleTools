//////////////////////////////////////////
// Global configuration
//////////////////////////////////////////

var CGI_URL = "http://ec2-79-125-52-174.eu-west-1.compute.amazonaws.com/cgi-bin/wiggleCGI.py?";
var attribute_values_file = "datasets.attribs.json";
var assembly = "GRCh37";

//////////////////////////////////////////
// Main function 
//////////////////////////////////////////

$(document).ready(
  function () {
    create_all_selectors();
    add_annotations();
    define_buttons();
  }
)

//////////////////////////////////////////
// Global variables
//////////////////////////////////////////

var selection_panels = [
  "choose",
  "chooseB",
  "chooseA2",
  "chooseA"
];

var panel_letters = {
  "choose":"A",
  "chooseB":"B",
  "chooseA2":"A",
  "chooseA":"A"
};

var attribute_values = null;

var reduction_opts = {
  "signal": {"Sum":"sum","Mininum":"min","Maximum":"max","Mean":"mean","Median":"median"},
  "regions": {"Intersection":"unit mult", "Union":"unit sum"}
};

var comparison_opts = {
  "regions": {"Intersection":"unit mult", "Union":"unit sum", "Difference": "unit diff"},
  "signal": {"Difference":"diff", "Ratio":"ratio", "Log Ratio":"log ratio", "T-test":"t-test", "Wilcoxon rank test":"wilcoxon"},
  "mixed": {"Distribution":"hist", "Profile curve":"profile", "Profile matrix":"profiles"}
};

//////////////////////////////////////////
// Creating multiselects 
//////////////////////////////////////////

function add_value_to_multiselect(value, div) {
  $("<option>").attr("value",value).text(value).appendTo(div);
}

function all_selects_are_used(panel) {
  var selects = panel.find(".form-control");
  var res = true; 
  selects.each(function(rank, select) {if ($(select).val() == "None") {res = false;}});
  return res;
}

function create_multiselect(container, attribute, panel) {
  var multiselect2 = $("<select>")
    .addClass("multiselect")
    .attr("multiple","multiple")
    .appendTo(container)
    .attr("attribute",panel_letters[panel.attr("id")]+ "_" + attribute);

  if (attribute in attribute_values) {
    attribute_values[attribute].map(function(value) {add_value_to_multiselect(value, multiselect2);});
  }
  multiselect2.multiselect({onChange: function(element, checked) {update_panel_count(panel);}, maxHeight: 400, buttonWidth:'100%'});
  multiselect2.parent().find('.btn').css("white-space","normal");
}

function reset_multiselect(select) {
  var panel = select.parents("[id*='choose']");
  var col = select.parent().parent().find(".multiselect").parent(".form-group");
  if (select.val() in attribute_values) {
    col.children().remove();
    create_multiselect(col, select.val(), panel);
  } else {
    col.parent().remove();
  }
  
  if (all_selects_are_used(panel)) {
    create_selection_div(panel);
  }
}

function create_attribute_select(container) {
  var select = $("<select>").addClass("form-control").appendTo(container);
  Object.keys(attribute_values).map(function(attribute) {add_attribute_to_select(attribute, select);});
  $("<option>").attr("value","None").text("None").attr("selected","selected").appendTo(select);
  select.change(function() {reset_multiselect(select);})
}

function add_attribute_to_select(attribute, select) {
  $("<option>").attr("value",attribute).text(attribute).appendTo(select);
}

function update_panel(panel) {
  // Compute initial count:
  update_panel_count(panel);
  // Set reduction select options
  update_panel_reduction(panel);
  // Comparison select:
  update_tab_comparison(panel.parents('.tab-pane'));
}

function create_selection_div(panel) {
  // Top most container
  var row = $("<div>").addClass("row").appendTo(panel.find("#selection"));

  // Division into fixed width columns
  var col1 = $("<div>").addClass("form-group").addClass("col-lg-5").appendTo(row);
  $("<div>").addClass("form-group").addClass("col-lg-1").text(" is: ").appendTo(row);
  var col2 = $("<div>").addClass("form-group").addClass("col-lg-5").appendTo(row);

  // Create attribute selector in column 1:
  create_attribute_select(col1);
  // Create empty value selector in column 2:
  create_multiselect(col2, "None", panel);

  update_panel(panel);
}

function create_all_selectors() {
  jQuery.getJSON(attribute_values_file).done(function(values) {
    attribute_values = values;
    selection_panels.map(function (id) {create_selection_div($("#"+id));});
  }).fail(catch_JSON_error);
}

//////////////////////////////////////////
// Creating the reference selectors 
//////////////////////////////////////////

function add_annotation(select, annotation) {
  $("<option>").text(name).appendTo(select);
}

function add_annotations_2(data) {
  var select = $('#refs');
  data['annotations'].map(function (annotation) {add_annotation(select, annotation)});
}

function add_annotations() {
  $.getJSON(CGI_URL + "assembly=" + assembly + "&annotations=1").done(add_annotations_2).fail(catch_JSON_error);
  fill_select($('#reference_reduction'), reduction_opts['regions']);
}

//////////////////////////////////////////
// Panel Query
//////////////////////////////////////////

function add_multi_select_query(list, multiselect) {
  var variable = multiselect.attr("attribute");
  var values = multiselect.val();  
  if (values != null && values != "") {
    values.map(function(value) {list.push([variable+'='+value])});
  }
  return list;
}

function  get_panel_type(panel) {
  return panel.find('#type').find('.active').find('input').attr('value');
}

function get_panel_letter(panel) {
  return panel_letters[panel.attr('id')];
}

function panel_query(panel) {
  list = []
  panel.find(".multiselect").each(function(i,value) {add_multi_select_query(list, $(value));});
  return list.join("&") + "&" + get_panel_letter(panel) + "_type=" + get_panel_type(panel);
}

//////////////////////////////////////////
// Computing panel selection count
//////////////////////////////////////////

function update_panel_count(panel) {
  url = CGI_URL + "count=true&assembly=" + assembly + "&" + panel_query(panel);
  $.getJSON(url).done(
   function(data, textStatus, jqXHR) {
     panel.find("#count").text(data["count"] + " elements selected");
   }
  ).fail(catch_JSON_error);
}

//////////////////////////////////////////
// Updating panel reduction select 
//////////////////////////////////////////

function fill_select(select, options) {
  select.children().remove();
  Object.keys(options).map(function(opt) {$("<option>").attr("value",options[opt]).text(opt).appendTo(select);});
}

function update_panel_reduction(panel) {
  fill_select(panel.find('#reduction'), reduction_opts[get_panel_type(panel)])
}

//////////////////////////////////////////
// Updating tab comparison select 
//////////////////////////////////////////

function update_tab_comparison(tab) {
  if (tab.attr('id') != 'compare' && tab.attr('id') != 'annotate') {
    return;
  } 
  var count = 0;  
  tab.find("[id*=choose]").each(function (index, panel) {if (get_panel_type($(panel)) == "regions") {count += 1;}});
  var select = tab.find("#comparison");

  // Sneaky way of detecting that we are in annotation tab, not comparison tab
  if (select.length == 0) {
    select = tab.find("#annotation");
    count += 1;
  }

  if (count == 2) {
    fill_select(select, comparison_opts["regions"]);
  } else if (count == 1) {
    fill_select(select, comparison_opts["mixed"]);
  } else {
    fill_select(select, comparison_opts["signal"]);
  }
}

//////////////////////////////////////////
// Button up!
//////////////////////////////////////////

function define_buttons() {
  $('#summary_button').click(summary);
  $('#comparison_button').click(comparison);
  $('#annotation_button').click(annotation);
  $('#result_button').click(get_result);
  $('[id*=type]').find('label').change(update_my_panel);
}

function update_my_panel() {
  $(this).addClass('active');
  var panel = $(this).parents('[id*="choose"]');
  var X = 1;
  update_panel(panel);
  $(this).removeClass('active');
}

function report_result(data) {
  if (data["status"] == "DONE") {
    var modal = $('#Success_modal').clone();
    modal.find('#url').attr('href',data['url']);
    modal.find('#view').attr('href',data['view']);
    modal.modal();
  } else if (data["status"] == "FAIL") {
    $('#Failure_modal').modal();	
  } else if (data["status"] == "UNKNOWN") {
    $('#Unknown_modal').modal();	
  } else {
    $('#Waiting_modal').modal();	
  }
}

// Get result
function get_result() {
  $.getJSON(CGI_URL + "result=" + $('#result_box').val()).done(report_result).fail(catch_JSON_error);
}

function return_ticket(data) {
  var modal = $("#JobSent_modal").clone();
  modal.find("#job_id").text(data["ID"]);
  modal.modal();
}

// Send job to server 
function submit_query(query) {
  $.getJSON(CGI_URL + "assembly=" + assembly + "&" + query).done(return_ticket).fail(catch_JSON_error);
}

// Request summary
function summary() {
  submit_query(panel_query($('#choose')) + '&wa=' + $('#choose').find('#reduction').val()); 
}

// Request comparison
function comparison() {
  submit_query(
    [ 
      panel_query($('#chooseA')),
      panel_query($('#chooseB')),
      'wa='+$('#chooseA').find('#reduction').val(),
      'wb='+$('#chooseB').find('#reduction').val(),
      'w='+$('#comparison').val()
    ].join("&")
  ); 
}

// Request annotation
function annotation() {
  submit_query(
    [ 
      panel_query($('#chooseA2')),
      $('#refs').val(),
      'wa='+$('#chooseA').find('#reduction').val(),
      'w='+$('#annotation').val()
    ].join("&")
  ); 
}

// JSON error handler
function catch_JSON_error(jqXHR, textStatus, errorThrown) {
  console.log('JSON failed: ' + textStatus + ":" + errorThrown + ":" + jqXHR.responseText + ":" + jqXHR.getAllResponseHeaders());
}
