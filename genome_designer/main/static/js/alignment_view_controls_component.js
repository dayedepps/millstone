/**
 * @fileoverview Component that decorates the controls for the view of a
 *     specific AlignmentGroup and the aligned ExperimentSamples.
 */

gd.AlignmentViewControlsComponent = gd.DataTableControlsComponent.extend({
  initialize: function() {
    gd.DataTableControlsComponent.prototype.initialize.call(this);

    this.render();
  },

  render: function() {
    this.decorateControls();
  },

  decorateControls: function() {
    this.drawDropdownOptions();
  },

  /** Draws dropdown options. */
  drawDropdownOptions: function() {
    // Option to delete samples.
    var downloadBamOptionHTML =
        '<a href="#" class="gd-id-ag-download-bam">Download BAM</a>';
    this.addDropdownOption(downloadBamOptionHTML);
    $('.gd-id-ag-download-bam').click(_.bind(this.handleDownloadBam, this));

    // Call SVs using contig assembly.
    if (FLAG_GENOME_FINISHING_ENABLED) {
      var assembleOptionHtml =
          '<a href="#" class="gd-id-assemble">Call SVs by de novo Assembly</a>';
      this.addDropdownOption(assembleOptionHtml);
      $('.gd-id-assemble').click(_.bind(this.handleAssembleContigs, this));
    }
  },

  handleDownloadBam: function() {
    var uidList = this.datatableComponent.getCheckedRowUids();

    // If nothing to do, show message.
    if (uidList.length != 1) {
      alert("Please select one sample at a time.");
      return;
    }

    var formJqueryObj = $('#gd-ag-download-bam-form');

    // Reset the form html
    formJqueryObj.empty();

    // Append the form fields.
    this._appendInputFieldToForm(formJqueryObj, 'estaUid', uidList[0]);

    // Submit the form. This cause a download to start.
    formJqueryObj.submit();
  },

  /** Helper method to append input value to form. */
  _appendInputFieldToForm: function(formJqueryObj, name, value) {
    formJqueryObj.append(_.template(
        '<input type="hidden" name="<%= name %>" value="<%= value %>">',
        {name: name, value: value}
    ));
  },

  /** Send request to generate contigs with default parameters **/
  handleAssembleContigs: function() {
    sample_alignment_uid_list = this.datatableComponent.getCheckedRowUids()

    var postData = {
        sampleAlignmentUidList: sample_alignment_uid_list
    };

    $.post('/_/alignments/generate_contigs', JSON.stringify(postData),
        _.bind(this.handleAssembleContigsResponse, this));
  },

  handleAssembleContigsResponse: function(response) {
    if (response.is_contig_file_empty == 1) {
      alert('No evidence for structural variants in this alignment');
    } else {
      this.trigger('MODELS_UPDATED');
    };
  },
});
